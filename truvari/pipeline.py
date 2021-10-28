"""
pipeline - Easy function chaining with futures parallelization
"""
import sys
import logging
import traceback

import concurrent.futures as cfuts
from functools import partial, reduce


def __reduce_wrapper(funcs, data, **kwargs):
    """
    This is a wrapper around reduce, which chains the functions
    We assume there must always be at least one argument
    """
    return reduce(lambda r, f: f(r, **kwargs), funcs, data)


def __partial_wrapper(funcs):
    """
    Wrap the reduce_wrapper in a partial so it may run multiple times
    (i.e. over multiple pieces of data)
    """
    return partial(__reduce_wrapper, funcs)


def fchain(pipe, data, workers=1):
    """
    Runs chained functions over every piece of data using a pool of futures.
    Yields results of chained functions called per data.

    pipe: list of [(function, kwargs), ...]
        - each function must be defined with at least `def foobar(data)`.
        - kwargs is a dictionary of other parameters the function may need.
        - if no parameters are needed, an item in the list can just be a function
    data: list data to process
        - data must be hashable
    workers: maximum number of ThreadPoolExecutor workers to create

    For example, the first yield of a 2 step pipe on the first piece of data is
    equivalent to:
        pipe[1][0](
            pipe[0][0](data[0], **pipe[0][1]),
            **pipe[1][1]
        )
    """
    m_reduce = __partial_wrapper([partial(_[0], **_[1]) if isinstance(_, tuple) else _ for _ in pipe])
    with cfuts.ThreadPoolExecutor(max_workers=workers) as executor:
        # Start the load operations and mark each future with its arg
        future_results = {executor.submit(m_reduce, _): _ for _ in data}
        for future in cfuts.as_completed(future_results):
            my_name = future_results[future]
            try:
                cur_result = future.result()
            except Exception as exc:  # pylint: disable=broad-except
                logging.critical('%r generated an exception: %s', my_name, exc)
                exc_type, exc_value, exc_traceback = sys.exc_info()
                logging.debug("Dumping Traceback:")
                traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stderr)
                # how do we stop safely... just return, or raise the exception?
                continue
            yield cur_result
            # Can't keep all of these around. Manually remove as you go.
            del(future_results[future])

def test():
    import pysam
    import time
    import random
    import itertools
    v = pysam.VariantFile("../repo_utils/test_files/input1.vcf.gz")
    def chunker(matcher, files):
        a = 40#random.randint(10, 50)
        for i in range(a):
            yield i, matcher, files
    def compare(null):
        # So this is the problem...
        number, matcher, mlist = null
        slp = random.randint(1, 3)
        time.sleep(slp)
        yield str(number) + '.' + str(slp) + '_' + str(matcher) + str(mlist)
    mat = "matcher"
    pipe = [compare]
    kwargs = {'matcher':mat, 'files':'balls'}#[('base', 1), ('comp', 2)] }
    c = chunker(mat, 'balls')
    for _ in itertools.chain.from_iterable(fchain(pipe, c, workers=8)):
        print(_)

if __name__ == '__main__':
    test()
