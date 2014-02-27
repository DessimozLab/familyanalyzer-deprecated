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
