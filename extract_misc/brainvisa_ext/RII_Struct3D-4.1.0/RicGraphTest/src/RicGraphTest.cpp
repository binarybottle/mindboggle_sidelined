// --------------------------- RicGraphTest.cpp -------------------------------
// Test for the RicGraph and RicGraphNode classes

#include <string>
#include "RicGraph.h"

using namespace std;

int main(int argc, char *argv[])
{

	// read a graph
	RicGraph TGraph("PH0040_Lgyri_default_session_auto.arg");

	cout << "Read in " << TGraph.nnodes << " nodes." << endl;

	// look for exact lablel name
	int n;
	if ( (n=TGraph.FindNodeNameExact("Pre-Central_left")) >= 0 )
	{
		cout << "Found " << "Pre-Central_left" << " at index " << n << endl;
	}
	else
	{
		cout << "Did not find Pre-Central_left" << endl;
	}

	// look for nodes with names containing a substring
	int cnt=0;
	cout << "Looking for Central" << endl;
	if ( (n=TGraph.FindNodeName("Central")) >= 0 )
	{
		cout << "Found " << TGraph.nodes[n].name << " at index " << n << endl;
		++cnt;
		while ( (n=TGraph.FindNextNodeName("Central")) >= 0 )
		{
			cout << "Found " << TGraph.nodes[n].name << " at index " << n << endl;
			++cnt;
		}
	}

	cout << "That's all folks!" << endl;

	return 0;
}

