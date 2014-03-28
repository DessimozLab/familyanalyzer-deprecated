try:
    from progressbar import ProgressBar, Percentage, Timer, ETA, Bar
    PROGRESSBAR = True
except ImportError:
    PROGRESSBAR = False

def setup_progressbar(msg, size):
    if not msg.endswith(': '):
        msg += ': '

    widgets = [msg,
               Percentage(), ' ',
               Bar(), ' ',
               Timer(), ' ',
               ETA()]

    pbar = ProgressBar(widgets=widgets, maxval=size)
    return pbar

def enum(*sequential, **named):
    """creates an Enum type with given values"""
    enums = dict(zip(sequential, range(len(sequential))), **named)
    enums['reverse'] = dict((value, key) for key, value in enums.items())
    return type('Enum', (object, ), enums)
