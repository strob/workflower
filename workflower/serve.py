# ui-based execution of a workflow
import json
import networkx as nx
import os
import subprocess
import traceback
from twisted.web.static import File
from twisted.web.resource import Resource
from twisted.web.server import Site
from twisted.internet import reactor

def wf2json(wf):
    nodes = {}                  # name -> subgraph or True
    edges = []

    for node in wf._graph.nodes():
        if hasattr(node, "_graph"):
            nodes[node.name] = wf2json(node)
        else:
            nodes[node.name] = False

    for edge in wf._graph.edges():
        edges.append((edge[0].name, edge[1].name))

    return {"nodes": nodes,
            "edges": edges}

def show_flow(wf):
    g = wf._graph
    for edge in g.edges():
        print edge, g.get_edge_data(*edge)["connect"]

def de_trait(t):
    out = {}
    for k,v in t.items():
        out[k] = str(getattr(t,k))
    return out

class NodeAPI(Resource):
    def __init__(self, node, nodeinfo):
        self.node = node
        self.nodeinfo = nodeinfo
        Resource.__init__(self)
    def render_GET(self, request):
        request.headers["Content-Type"] = "application/json"
        return json.dumps([X.get_obj() for X in self.nodeinfo])

class ExecutionInfo:
    def set_inputs(self, inputs):
        self.inputs = inputs
    def set_outputs(self, outputs):
        self.outputs = outputs
    def get_obj(self):
        out = {}
        if hasattr(self, "inputs"):
            out["inputs"] = de_trait(self.inputs)
        if hasattr(self, "outputs"):
            out["outputs"] = de_trait(self.outputs)
        return out

class WorkflowAPI(Resource):
    def __init__(self, wf):
        self.wf = wf
        self.nodeinfo = {}      # node_name -> [ExecutionInfo]
        Resource.__init__(self)

    def getChild(self, name, req):
        return NodeAPI(self.wf.get_node(name), self.nodeinfo.get(name, []))

    def run(self, graph, *a, **kw):
        "Custom workflow plugin to track execution"
        print "graph", graph, graph.nodes()
        for node in nx.topological_sort(graph):
            print "running!", node, node.name
            ei = ExecutionInfo()
            ei.set_inputs(node.inputs)
            self.nodeinfo.setdefault(node.name, []).append(ei)

            try:
                res = node.run()
            except Exception, err:
                print err
                self.runout = {"type": "error",
                           "node": str(node.name),
                           "error": traceback.format_exc()}
                return
            ei.set_outputs(res.outputs)
            self.runout = {"type": "success"}

    def render_GET(self, request):
        query = request.args.get("q", [])
        if len(query) > 0:
            if query[0] == 'json':
                request.headers["Content-Type"] = "application/json"
                return json.dumps(wf2json(self.wf))

            elif query[0] == 'dot':
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
                        colored=False)))
                return stdout
            elif query[0] == "run":
                request.headers["Content-Type"] = "application/json"
                self.wf.run(plugin=self)
                return json.dumps(self.runout)
            elif query[0] == "shell":
                import pdb; pdb.set_trace()
                return "hope you fixed it!"
                
        return "The workflow is called: %s" % (self.wf.name)

def run(wf, port=9090):
    # web-control of a workflow
    datapath = os.path.join(os.path.dirname(__file__), "data")
    root = File(datapath)
    workflow_api = WorkflowAPI(wf)
    root.putChild("wf", workflow_api)
    site = Site(root)

    reactor.listenTCP(port, site, interface="0.0.0.0")
    print 'listening on port '+str(port)
    reactor.run()
    

if __name__=='__main__':
    import sample_workflow
    # show_flow(sample_workflow.wf)
    # sample_workflow.wf.run(plugin=MRun())

    run(sample_workflow.wf)
