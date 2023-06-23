"""
\file Model.cpp \brief Storage of a model, its restraints,
                        constraints and particles.

Copyright 2007-2022 IMP Inventors. All rights reserved.
"""

class Model:
    def __init__(self, name="Model1"):
        self.name = name
        self.moved_particles_restraint_cache_ = None
        self.moved_particles_particle_cache_ = None
        self.cur_stage_ = 0 # NOT_EVALUATING
        self.first_call_ = True
        self.age_counter_ = 1
        self.dependencies_age_ = 0
        self.saved_dependencies_age_ = 0
        self.dependencies_saved_ = False
        self.moved_particles_cache_age_ = 0

    class Masks:
        read_mask_ = []
        write_mask_ = []
        add_remove_mask_ = []
        read_derivatives_mask_ = []
        write_derivatives_mask_ = []

    class ScoreStates:
        def __init__(self):
            self.score_states = []

        def append(self, score_state):
            self.score_states.append(score_state)
