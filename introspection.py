import importlib
import os

def module_walk(root, package=None):
    "yield all child-modules of a parent"

    yield root
    siblings = os.listdir(os.path.dirname(root.__file__))
    if package is None:
        package = root.__package__

    for script in filter(lambda x: x.endswith(".py"), siblings):
        if not script.startswith("__"):
            modname = "%s.%s" % (package, script[:-3])
            yield importlib.import_module(modname, package=package)

    for subdir in filter(lambda x: os.path.isdir(
            os.path.join(
                os.path.dirname(root.__file__), x)),
                         siblings):

        modname = "%s.%s" % (package, subdir)
        try:
            parent_module = importlib.import_module(modname, package=package)
        except ImportError:
            # XXX: What about subdirectories lacking an __init__.py?
            continue
        yield parent_module
        for m in module_walk(parent_module, package=modname):
            yield m

def find_subclasses(module, parent_class):
    yielded = set()             # eliminate probable duplicate
    for m in module_walk(module):
        for c in dir(m):
            obj = getattr(m, c)
            if not c in yielded and type(obj) == type and issubclass(obj, parent_class):
                yield obj
                yielded.add(c)

if __name__=='__main__':
    import nipype
    import nipype.interfaces.base
    import nipype.pipeline.engine

    print "INTERFACES"

    for x in find_subclasses(nipype, nipype.interfaces.base.Interface):
        print x.__name__

    print "WORKFLOWS"

    for y in find_subclasses(nipype, nipype.pipeline.engine.WorkflowBase):
        print y.__name__

