#include <fstream>
#include <string>
#include <iostream>
#include <vector>
using namespace std;

int nNode;
int nLink;
int nReq;
int nPath;

vector<int> linkCapa;			// capacity of all links, including reverse links
vector<int> reqBw;				// bandwidth demand of all requests
vector<vector<vector<int>>> reqPaths;			// alternative paths for all requests

// ×Ö·û´®·Ö¸îº¯Êý  
vector<string> split(string str, string pattern)
{
	string::size_type pos;
	vector<string> result;

	str += pattern;
	int size = str.size();

	for (int i = 0; i < size;)
	{
		pos = str.find(pattern, i);
		if (pos == i)									// avoid space in the end
			break;
		string s = str.substr(i, pos - i);
		result.push_back(s);
		i = pos + pattern.size();
	}
	return result;
}

// read topo and requests data

void readTopoAndRequest(const char *filename)
{
	ifstream inputfile(filename);
	vector<string> strs;
	string s;
	
	// number of nodes and links
	getline(inputfile, s);				
	strs = split(s, " ");
	nNode = stoi(strs[0]);
	nLink = stoi(strs[1]);

	// capacity of all links, including reverse links
	linkCapa.reserve(nLink * 2);
	for (int i = 0; i < nLink * 2; ++i)
		linkCapa.push_back(0);

	for (int i = 0; i < nLink; ++i)
	{
		getline(inputfile, s);
		strs = split(s, " ");
		int start = stoi(strs[0]);
		int end = stoi(strs[1]);
		int capa = stoi(strs[2]);
		if (end == start + 1)		// vertical link
		{
			linkCapa[start / 20 * 19 + start % 20] = capa;		// arrange the order of links: first vertical links, then horizontal links
			linkCapa[nLink + start / 20 * 19 + start % 20] = capa;          // reverse link
		}
		else
		{
			linkCapa[19 * 25 + start % 20 * 24 + start / 20] = capa;
			linkCapa[nLink + 19 * 25 + start % 20 * 24 + start / 20] = capa;
		}
	}

	// number of requests and redundant paths
	getline(inputfile, s);
	strs = split(s, " ");
	nReq = stoi(strs[0]);
	nPath = stoi(strs[1]);

	// requests and paths
	reqBw.reserve(nReq);
	reqPaths.reserve(nReq);

	for (int i = 0; i < nReq; ++i)
	{
		// bandwidth request
		getline(inputfile, s);
		strs = split(s, " ");
		reqBw.push_back(stoi(strs[1]));

		// all paths
		vector<vector<int>> paths;
		for (int j = 0; j < nPath; ++j)
		{
			vector<int> path;
			getline(inputfile, s);
			strs = split(s, " ");

			int preNode = stoi(strs[0]);
			for (int k = 1; k < strs.size(); ++k)
			{
				int curNode = stoi(strs[k]);
				if (curNode == preNode + 1)		// from up to down
					path.push_back(preNode / 20 * 19 + preNode % 20);        // append the link index
				else if (curNode == preNode - 1)            // from down to up
					path.push_back(nLink + preNode / 20 * 19 + preNode % 20);
				else if (curNode == preNode + 20)           // from left to right
					path.push_back(19 * 25 + preNode % 20 * 24 + preNode / 20);
				else                              // from right to left
					path.push_back(nLink + 19 * 25 + preNode % 20 * 24 + preNode / 20);
				preNode = curNode;
			}
			paths.push_back(path);
		}
		reqPaths.push_back(paths);
	}

	nLink *= 2;			// include reverse paths

}


