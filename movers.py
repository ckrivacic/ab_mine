from pyrosetta import *
from mover_utils import *


def get_minimizer(residues_bb_movable, residues_sc_movable):
    mm = rosetta.core.kinematics.MoveMap()
    mm.set_chi(False)
    mm.set_bb(False)

    for i in residues_bb_movable:
        mm.set_bb(i, True)
        mm.set_chi(i, True)

    for i in residues_sc_movable:
        mm.set_chi(i, True)

    min_opts = rosetta.core.optimization.MinimizerOptions(
            "lbfgs_armijo_nonmonotone", 0.01, True )


    minmover = rosetta.protocols.minimization_packing.MinMover()
    minmover.movemap(mm)
    minmover.min_options(min_opts)
    sfxn = create_score_function('ref2015_cst')
    minmover.score_function(sfxn)

    return minmover


def get_loop_modeler(pose, residues, movemap=None, task_factory=None, 
        mover='ngk', fast=False, resbuffer=3):
    '''Run loop modeler on the pose (default to NGK)
    Available movers: 
        - NGK
        - Loophash KIC'''
    '''
    Note for testing loop modeling parameters:
    To change temp/number of cycles/etc., get the centroid stage and
    fullatom stage LoopProtocol objects via loopmodeler.centroid_stage()
    and loopmodeler.fullatom_stage(). Can set temp and cycles from
    there.
    '''
    # Get start and end of loop from list of residues
    start_residue = residues[0]
    end_residue = residues[-1]
    # Simple setup movemap
    mm = setup_movemap(residues, residues)
    sfxn = setup_restrained_sfxn(['coordinate_constraint'],[1.0])

    loopmodeler = rosetta.protocols.loop_modeler.LoopModeler()
    if mover=='ngk':
        loopmodeler.setup_kic_config()
        #loops = generate_loops_from_res_selector(pose, designable_selector,
        #        focus_residue, resbuffer=resbuffer)
        loops = generate_loops_simple(pose, start_residue, end_residue)
    elif mover=='lhk':
        '''A note on LHK: You need to mutate focus residues to their motif
        residue before running.'''
        #assert(resbuffer >= 4)
        loopmodeler.setup_loophash_kic_config(True, '')
        loops = generate_loops_simple(pose, start_residue, end_residue)


    loopmodeler.set_loops(loops)
    loopmodeler.set_fa_scorefxn(sfxn)

    if fast:
        loopmodeler.centroid_stage().mark_as_test_run()
        loopmodeler.fullatom_stage().mark_as_test_run()

    if task_factory:
        loopmodeler.set_task_factory(task_factory)

    return loopmodeler
