import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe


ii1 = util.IdentityInterface(fields=['a'])
node_1 = pe.Node(ii1, name="ii_1")
node_1.iterables = ("a", [42,24])

def double(a):
    return 2*a

fn1 = util.Function(input_names=['b'],
                    output_names=['c'],
                    function=double)
node_2 = pe.Node(fn1, name="fn1")

wf = pe.Workflow(name="sample_workflow")
wf.connect(node_1, 'a', node_2, 'b')
# wf.run()
