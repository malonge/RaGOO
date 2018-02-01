import time
import subprocess


def run(cmnd):
    """ Run command and report status. """
    log(' ---- Running : %s' % cmnd)
    if subprocess.call(cmnd, shell=True) != 0:
        raise RuntimeError('Failed : %s ' % cmnd)


def log(message):
    """ Log messages to standard output. """
    print (time.ctime() + ' ' + message)