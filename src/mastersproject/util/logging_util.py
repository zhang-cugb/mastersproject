# --- LOGGING UTIL ---
import logging
from decorator import decorator
import time


@decorator
def warn_slow(func, timelimit=60, *args, **kw):
    t0 = time.time()
    result = func(*args, **kw)
    dt = time.time() - t0
    if dt > timelimit:
        logging.warning('%s took %d seconds', func.__name__, dt)
    else:
        logging.info('%s took %d seconds', func.__name__, dt)
    return result


@decorator
def timer(func, *args, **kw):
    t0 = time.time()
    result = func(*args, **kw)
    dt = time.time() - t0
    logging.info(f'{func.__name__} took {dt} seconds')
    return result


@decorator
def trace(f, *args, **kw):
    kwstr = ', '.join('%r: %r' % (k, kw[k]) for k in sorted(kw))
    logging.debug(f"Calling {f.__name__} with args {args}, {{{kwstr}}}")
    return f(*args, **kw)