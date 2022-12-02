from functools import wraps
from time import time


def timing(f):
    """Decorator that allows the timing of function using the wrapping functool."""

    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print(f"func:{f.__name__} args:[{args}, {kw}] took: {te-ts:2.4f} sec")
        return result

    return wrap
