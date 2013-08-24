# ui-based execution of a workflow
import json
import networkx as nx
import subprocess
from twisted.web.static import File
from twisted.web.resource import Resource
from twisted.web.server import Site
from twisted.internet import reactor

def show_flow(wf):
    g = wf._graph
    for edge in g.edges():
        print edge, g.get_edge_data(*edge)["connect"]

class MRun:
    def run(self, graph, *a, **kw):
        print "my_runner"
        print "graph", graph, graph.nodes()
        for node in nx.topological_sort(graph):
            print "running!", node
            res = node.run()
            print "result", res, res.outputs

def de_trait(t):
    out = {}
    for k,v in t.items():
        out[k] = str(getattr(t,k))
    return out

class NodeAPI(Resource):
    def __init__(self, node):
        self.node = node
        Resource.__init__(self)
    def render_GET(self, request):
        request.headers["Content-Type"] = "application/json"
        return json.dumps({"inputs": de_trait(self.node.inputs),
                           "outputs": de_trait(self.node.outputs)})

class WorkflowAPI(Resource):
    def __init__(self, wf):
        self.wf = wf
        Resource.__init__(self)

    def getChild(self, name, req):
        return NodeAPI(self.wf.get_node(name))

    def render_GET(self, request):
        query = request.args.get("q", [])
        if len(query) > 0:
            if query[0] == 'dot':
                request.headers["Content-Type"] = "text/plain"
                return self.wf._get_dot()
            elif query[0] == 'svg':
                request.headers["Content-Type"] = "application/svg"

                p = subprocess.Popen(["dot", "-Tsvg"], 
                                     stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE)
                stdout, stderr = p.communicate(
                    input="digraph %s{\n%s}" % (self.wf.name, self.wf._get_dot(
                        simple_form=True,
                        colored=True)))
                return stdout
            elif query[0] == "run":
                self.res = self.wf.run()
                return "OK!"
            elif query[0] == "shell":
                import pdb; pdb.set_trace()
                return "hope you fixed it!"
                
        return "The workflow is called: %s" % (self.wf.name)

def run(wf):
    # web-control of a workflow
    root = File("www")
    workflow_api = WorkflowAPI(wf)
    root.putChild("wf", workflow_api)
    site = Site(root)

    reactor.listenTCP(9090, site, interface="0.0.0.0")
    print 'listening on port 9090'
    reactor.run()
    

if __name__=='__main__':
    import sample_workflow
    # show_flow(sample_workflow.wf)
    # sample_workflow.wf.run(plugin=MRun())

    run(sample_workflow.wf)
