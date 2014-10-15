import os

def data_dir():
    """
    Returns the data directory that contains files for data and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'data')


def test_dir():
    
    """
    Returns the data directory that contains example files for tests and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'test', 'data')

def data_file(fn):
    fn = os.path.join(data_dir(), fn)

    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % (fn))
    return fn

def test_file(fn):
    fn = os.path.join(test_dir(), fn)
    
    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % (fn))
    return fn
