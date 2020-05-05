#! /wynton/home/kortemme/krivacic/software/anaconda36/bin/python3
#$ -l mem_free=4G
#$ -cwd

from pyrosetta import *
from movers import *
from mover_utils import *
from utils import *
import glob
import pickle as pkl
import pandas as pd

def score_ab(pose, minimize=False, kic=False, random=False, lhk=False):
    resnums = [95,96,97,98,99,100]
    pdbinfo = pose.pdb_info()
    rosettastart = pdbinfo.pdb2pose('H', resnums[0])
    chain = pose.chain(rosettastart)
    chainpose = pose.split_by_chain(chain)
    sfxn = create_score_function('ref2015')
    rosetta_resnums = rosetta_numbers_from_pdb(resnums, chainpose,
            chain='H')
    rosetta_resnums = [x for x in rosetta_resnums if x!=0]
    if len(rosetta_resnums) < 6:
        return 0

    # Setup fold tree
    ft = FoldTree()
    ft.add_edge(1, rosetta_resnums[0], -1)
    middle = (rosetta_resnums[-1] + rosetta_resnums[0])//2
    ft.add_edge(rosetta_resnums[0], middle, -1)
    ft.add_edge(rosetta_resnums[0], rosetta_resnums[-1], 1)
    ft.add_edge(rosetta_resnums[-1], middle + 1, -1)
    ft.add_edge(rosetta_resnums[-1], chainpose.size(), -1)

    print('Setting up the following fold tree:')
    print(ft)
    ft.check_fold_tree()
    chainpose.fold_tree(ft)

    if kic:
        loopmodeler = get_loop_modeler(chainpose, rosetta_resnums,
                fast=True)
        loopmodeler.apply(chainpose)

    if lhk:
        loopmodeler = get_loop_modeler(chainpose, rosetta_resnums,
                mover='lhk', fast=True)
        loopmodeler.apply(chainpose)

    if random:
        #loopprotocol = rosetta.protocols.loop_modeling.LoopProtocol()
        kicmover = rosetta.protocols.kinematic_closure.KicMover()
        #loopprotocol.add_mover(kicmover)
        start_residue = rosetta_resnums[0]
        end_residue = rosetta_resnums[-1]
        loop = generate_loops_simple(chainpose, start_residue,
                end_residue)
        kicmover.set_loops(loop)
        kicmover.apply(chainpose)


    if minimize:
        minimizer = get_minimizer(rosetta_resnums, rosetta_resnums)
        minimizer.apply(chainpose)

    total = sfxn(chainpose)

    loop_energy = 0
    residue_dict = {}
    ppo_vector = []
    phipsi_vector = []
    for resnum in resnums:
        rosettanum = chainpose.pdb_info().pdb2pose('H', resnum)
        ppo_vector.append(chainpose.phi(rosettanum))
        ppo_vector.append(chainpose.psi(rosettanum))
        ppo_vector.append(chainpose.omega(rosettanum))
        phipsi_vector.append(chainpose.phi(rosettanum))
        phipsi_vector.append(chainpose.psi(rosettanum))

        if rosettanum != 0:
            residue_energy=\
                    chainpose.energies().residue_total_energy(rosettanum)
            loop_energy += residue_energy
            residue_dict[rosettanum] = residue_energy
        else:
            print('Error: PDB number {}, chain H, has no associated'\
                        ' rosetta number.'.format(resnum))
    
    return {'pdbid': pdbinfo.name(), 
            'minimized': minimize,
            'kic': kic,
            'lhk': lhk,
            'random': random,
            'res_scores': [residue_dict],
            'pose_score': total,
            'loop_score': loop_energy,
            'ppo': ppo_vector,
            'phipsi': phipsi_vector}

def test_score():
    init('-ignore_unrecognized_res -ex1 -ex2 ' +\
            '-lh:db_path=/wynton/home/kortemme/krivacic/rosetta/database/loophash_db '+\
            '-lh:loopsizes 6 ' + \
            '-total_threads 1')
    pose = pose_from_pdb('pdbs/1BZQ_1.pdb')
    out0 = score_ab(pose, lhk=True)
    print(out0)
    out1 = score_ab(pose, minimize=False)
    out2 = score_ab(pose, minimize=True)
    out3 = score_ab(pose, kic=True)
    out4 = score_ab(pose, random=True)
    print(out1)
    print(out2)
    print(out3)
    print(out4)


def cluster_run():
    tasknum = int(os.environ['SGE_TASK_ID']) - 1
    #tasknum = 621
    pdbs = sorted(glob.glob('pdbs/*.pdb'))
    # 50 pdbs per task
    start = tasknum * 5
    end = tasknum * 5 + 5

    rows = []
    for i in range(start, end):
        if i < len(pdbs):
            pdb = pdbs[i]
        else: break
        print('Running pdb {}'.format(pdb))
        for protocol in ['minimize', 'score', 'kic', 'lhk', 'random']:
        
            minimize = False
            kic = False
            lhk = False
            random = False

            if protocol=='score':
                num_iter = 1
            elif protocol=='minimize':
                minimize = True
                num_iter = 1
            elif protocol=='kic':
                kic = True
                num_iter = 10
            elif protocol=='lhk':
                lhk = True
                num_iter = 20
            elif protocol=='random':
                random = True
                num_iter = 30

            for j in range(0, num_iter):
                pose = pose_from_file(pdb)
                print('Running protocol {} on {}, iter {}'.format(protocol, 
                        pdb, j))
                #try:
                row = score_ab(pose, minimize=minimize, 
                        kic=kic, random=random, lhk=lhk)
                print(row)
                if row==0:
                    print('PDB {} had fewer than 6 residues in loop'.format(pdb))
                    continue
                else:
                    rows.append(row)
                #except:
                #    print('ERROR running {} using protocol {}'.format(pdb, protocol))
    
    df = pd.DataFrame(rows)
    df.to_pickle('out/df_{}.pkl'.format(tasknum))

if __name__=='__main__':
    #test_score()
    init('-ignore_unrecognized_res -ex1 -ex2 ' +\
            '-lh:db_path=/wynton/home/kortemme/krivacic/rosetta/database/loophash_db '+\
            '-lh:loopsizes 6 ' + \
            '-total_threads 1')
    cluster_run()
