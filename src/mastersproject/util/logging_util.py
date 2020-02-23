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


def __setup_logging(path, log_fname="results.log"):
    path = str(path)
    # GTS logger
    gts_logger = logging.getLogger('GTS')
    gts_logger.setLevel(logging.INFO)
    if gts_logger.handlers:
        gts_logger.handlers = []

    # PorePy logger
    pp_logger = logging.getLogger('porepy')
    pp_logger.setLevel(logging.DEBUG)
    if pp_logger.handlers:
        pp_logger.handlers = []

    common_format = logging.Formatter(logging.BASIC_FORMAT)

    # Add handler for logging debug messages to file.
    fh = logging.FileHandler(path + "/" + log_fname, 'w+')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(common_format)

    gts_logger.addHandler(fh)
    pp_logger.addHandler(fh)

    # Add handler for logging info messages to console
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(common_format)

    gts_logger.addHandler(ch)
    pp_logger.addHandler(ch)




