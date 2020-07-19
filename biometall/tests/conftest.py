import os

TESTPATH = os.path.dirname(os.path.abspath(__file__))

def datapath(path):
    return os.path.join(TESTPATH, 'data', path)