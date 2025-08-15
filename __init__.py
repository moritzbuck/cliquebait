
_verbose_level_=1
version = 0.1
default_gap_size = 2.5
default_strain_cutoff = 99.9
default_lower_cutoff = 89.0
default_denoising_cutoff = 0.95

def set_verbose(level):
    _verbose_level_ = level

def get_verbose():
    return _verbose_level_ 
