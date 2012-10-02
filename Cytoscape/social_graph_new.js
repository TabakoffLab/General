// Visualize Relationship in Twitter by using Google Social Graph API

importPackage( Packages.javax.swing );
importPackage(java.io);
importPackage(java.net);
importPackage( Packages.cytoscape.layout );
importPackage( Packages.cytoscape );

// Base URL of Twitter
var twitterURL = "http://twitter.com/";

// Change this to your id
var myID = "c_z";

// Depth of search
// As you know, in this type of social networks, 
// number of nodes grows exponentially if you increase this number!!
var depth = 1;


var nodeAttr = Cytoscape.getNodeAttributes();
var newNetwork = Cytoscape.createNetwork("Twitter Graph from " + myID);

function expand(id) {
	var me = Cytoscape.getCyNode(id, true);
	newNetwork.addNode(me);

	var nodes = [];
	var url = new URL('http://socialgraph.apis.google.com/lookup?q=' + twitterURL + id + '&edo=1&edi=1');
	var line, json = '';

	try {
		var stream = new BufferedReader(new InputStreamReader(url.openStream()));

		while(line = stream.readLine()) 
			json += line;  

		stream.close();
	}catch (ex) {
		var empty = [];
		return empty;
	}
	
	relations = eval("(" + json + ")");

	var people = relations.nodes[twitterURL+id].nodes_referenced;
	
	var nodeIDs = [];
	var idx = 0;
	for(key in people) {
		if(people[key]['types'] != 'me') {
			var newid = key.substring(19, key.length);
			var targetNode = Cytoscape.getCyNode(newid, true);
			nodeAttr.setAttribute(newid, "Twitter URL", twitterURL+newid);
			newNetwork.addNode(targetNode);
			nodes.push(targetNode);
			nodeIDs[idx] = newid;
			idx++;
		}
	}
	
	for(i=0; i<nodes.length; i++) {
		edge = Cytoscape.getCyEdge(me, nodes[i], "interaction", "is_following", true);
		if(edge != null)
			newNetwork.addEdge(edge);
	}
	
	var followers = relations.nodes[twitterURL+id].nodes_referenced_by;
	
	nodes = [];
	idx = 0;
	for(key in followers) {
		if(followers[key]['types'] != 'me') {
			var newid = key.substring(19, key.length);
			var targetNode = Cytoscape.getCyNode(newid, true);
			nodeAttr.setAttribute(newid, "Twitter URL", twitterURL+newid);
			newNetwork.addNode(targetNode);
			nodes.push(targetNode);
			nodeIDs.concat(newid);
			idx++;
		}
	}
	
	for(i=0; i<nodes.length; i++) {
		edge = Cytoscape.getCyEdge(me, nodes[i], "interaction", "is_followed_by", true);
		if(edge != null)
			newNetwork.addEdge(edge);
	}
	return nodeIDs;
}
	
var idBuffer = new Array();
var buffers = [];
buffers.push(expand(myID));

for(counter=0; counter<depth; counter++) {

	var buffers2 = [];
	for(k=0; k<buffers.length; k++) {
		var idBuffer = buffers.pop(); 
		for(j=0; j<idBuffer.length; j++) {
			var nextID = idBuffer[j];
			if(nextID != null)
				buffers2.push(expand(nextID));
		}
	}
	for(l=0; l<buffers2.length; l++)
		buffers.push(buffers2.pop());
}