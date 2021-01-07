/*----------------------------------------------------------------*/
/* INCLUDES */
/*----------------------------------------------------------------*/

//OBJ Loader
#include "OBJ_Loader.h"

using namespace std;
using namespace objl;

//Custom structs
#include "structs.h"
using namespace structs;

//Library for solving cubic equations
//Library can be found here http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
//#include "poly34.cpp"

//standard includes
#include <math.h>
#include <sys/stat.h>
#include <string>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <omp.h>
#include <chrono> 
#include <string.h>
using namespace std::chrono; 


/*----------------------------------------------------------------*/
/* DEFINITIONS */
/*----------------------------------------------------------------*/

#define HELP    "Usage: %s -f <filePath> -n <nodeCount> [-ci]"\
                "\n"\
                "\n -f <filePath>: .obj input file"\
                "\n -n <nodeCount>: number of nodes to generate"\
                "\n -c: honeycomb mode"\
                "\n -i: improve node grid (remove nodes close to object surface)"\
                "\n"\
                "\n Generates a node grid inside an Axis-Aligned Bounding Box of an object and exports it into a .node file."\
                "\n Grid properties:"\
                "\n 	- size of the  grid is adjusted to object size"\
                "\n 	- nodes can be outside of the object"\
                "\n 	- grid is aligned with the object x,y,z axes"\
                "\n"\
                "\n"

/*----------------------------------------------------------------*/
/* GLOBAL VARS */
/*----------------------------------------------------------------*/


/*----------------------------------------------------------------*/
/* FUNCTIONS */
/*----------------------------------------------------------------*/

void outputDuration(double dcount)
{
	if(dcount > 1000000)
	{
		dcount /= 1000000;
		cout << " " << dcount << " s"; 
	}
	else if(dcount > 1000)
	{
		dcount /= 1000;
		cout << " " << dcount << " ms"; 
	}
	else
	{
		cout << " " << dcount << " \xC2\xB5s"; 
	}
}

/******************************************************************
* Create a 3D cubic node grid and save it to a .node file
* Input: 	center 		center of the grid
			length		edge length of the cubic grid (x,y,z)
			count 		node count to generate
			outputDir	output directory for node file
			fileName	name of the xxx.node file
* Output:	xxx.node file
* Return:
*******************************************************************/
void generateNodeGrid(OBJ &obj, int count, string outputDir, string fileName, bool honeycombMode, bool improve)
{
	/******************************************************************
	**** STEP SIZE CALCULATION
	******************************************************************/
	//execution time measuring
	auto timeStart = high_resolution_clock::now();
	auto timeStop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(timeStop - timeStart); 

	//step sizes
	double averageStepSize;
	double stepX, stepY, stepZ;
	

	//calculate approximate step size
	averageStepSize = cbrt(obj.edgeX*obj.edgeY*obj.edgeZ)/(cbrt(count) + 1);
	cout << "=> Approximate Step Size: " << averageStepSize << endl;

	//calculate nodes count per dimension
	int nodesX = round(obj.edgeX / averageStepSize) - 1;
	int nodesY = round(obj.edgeY / averageStepSize) - 1;
	int nodesZ = round(obj.edgeZ / averageStepSize) - 1;

	cout << "Exact Node Count (x,y,z): " << obj.edgeX / averageStepSize << ", " << obj.edgeY / averageStepSize << ", " << obj.edgeZ / averageStepSize << endl;

	//at least one node per dimension (x,y,z)
	if(nodesX < 1) nodesX = 1;
	if(nodesY < 1) nodesY = 1;
	if(nodesZ < 1) nodesZ = 1;

	//calculate step size per dimension
	stepX = obj.edgeX / (nodesX + 1);
	stepY = obj.edgeY / (nodesY + 1);
	stepZ = obj.edgeZ / (nodesZ + 1);

	averageStepSize = (stepX + stepY + stepZ)/3.0;


	cout << endl;

	cout << " --- Circumscribing Rectangle Cuboid --- " << endl;
	cout << "Upper Boundary Corner (x,y,z): " << obj.maxX << ", " << obj.maxY << ", " << obj.maxZ << endl;
	cout << "Lower Boundary Corner (x,y,z): " << obj.minX << ", " << obj.minY << ", " << obj.minZ << endl;
	cout << "Edge Lengths (x,y,z): " << obj.edgeX << ", " << obj.edgeY << ", " << obj.edgeZ << endl;

	cout << endl << " --- Node Grid Parameters --- " << endl;
	cout << "Average Step Size: " << averageStepSize << endl;
	cout << "Step Sizes (x,y,z): " << stepX << ", " << stepY << ", " << stepZ << endl;
	cout << "Nodes Count (x,y,z): " << nodesX << ", " << nodesY << ", " << nodesZ << endl;
	cout << "Honeycomb Mode: " << honeycombMode << endl;
	cout << "Mesh Triangle Count: " << obj.triangleCount << endl;

	cout << endl;

	/******************************************************************
	**** NODE GRID GENERATION
	******************************************************************/
	//honeycomb mode
	int honeycombToggle = 0;
	double honeycombOffsetY = stepY/2.0;
	double honeycombOffsetZ = stepZ/2.0;

	//calculate first node position
	double gridX = obj.minX + stepX;
	double gridY = obj.minY + stepY;
	double gridZ = obj.minZ + stepZ;

	//initialize double linked list with first grid node
	Node *curNode = new Node();

	//generate Node Grid
	unsigned int nodeCount = 0;

	int nodeCountX = 0;
	int nodeCountY = 0;
	int nodeCountZ = 0;

	cout << endl << "Generating Node Grid..." << endl;

	timeStart = high_resolution_clock::now();
	//x loop
	while(nodeCountX < nodesX)
	{
		//y loop
		while(nodeCountY < nodesY)
		{
			//z loop
			while(nodeCountZ < nodesZ)
			{
				curNode = curNode->appendNode(new Node(Vector(gridX, gridY, gridZ)));

				nodeCountZ++;
				gridZ += stepZ;

				nodeCount++;
			}

			nodeCountY++;
			gridY += stepY;

			nodeCountZ = -honeycombToggle;
			gridZ = obj.minZ + stepZ - honeycombOffsetZ*honeycombToggle;
		}

		//honyecomb mode toggle
		if(honeycombMode) honeycombToggle = 1 ^ honeycombToggle;

		nodeCountX++;
		gridX += stepX;

		nodeCountY = -honeycombToggle;
		gridY = obj.minY + stepY - honeycombOffsetY*honeycombToggle;
		nodeCountZ = -honeycombToggle;
		gridZ = obj.minZ + stepZ - honeycombOffsetZ*honeycombToggle;
	}

	timeStop = high_resolution_clock::now();
	duration = duration_cast<microseconds>(timeStop - timeStart); 
  
	cout << "Grid Node Count: " << nodeCount << endl;

	cout << "Node Grid Generation Time: ";
    outputDuration(duration.count());
    cout << endl;

    /******************************************************************
	**** NODE GRID IMPROVEMENT
	******************************************************************/
	if(improve)
	{
		//remove nodes too close to triangle mesh
		//take smallest step size and factor in double precision error
		double nodeRemoveThreshold = (stepX < stepY) ? stepX : stepY;
		nodeRemoveThreshold = (stepZ < nodeRemoveThreshold) ? stepZ : nodeRemoveThreshold;
		nodeRemoveThreshold = nodeRemoveThreshold/2.0 - nodeRemoveThreshold/1000.0;

		cout << "Removing Nodes close to surface mesh..." << endl;
		cout << "Node Removal Threshold Distance: " << nodeRemoveThreshold << endl;

		//make array with pointers to double-linked for easier parallelization
		Node *nodeArray[nodeCount];
		curNode = curNode->getFirst();

		for(unsigned int i = 0; i < nodeCount; i++)
		{
			nodeArray[i] = curNode;
			if(curNode->nextNode != NULL) curNode = curNode->nextNode;
		}

		//Parallel Loop marking nodes for deletion
		timeStart = high_resolution_clock::now();
		#pragma omp parallel for
		for(unsigned int i = 0; i < nodeCount; i++)
		{
			for(unsigned int j = 0; j < obj.triangleCount; j++)
			{
				Vector closestPoint = obj.triangles[j].ClosestPoint(nodeArray[i]->p);

				//Sanity Check for point being outside of triangle - deactivated by default for better performance
				//if(!obj.triangles[i].CheckCollision(closestPoint)) cout << "Point Outside Triangle?!" << endl;
				
				double distanceToTriangle = (nodeArray[i]->p - closestPoint).Length();

				if(distanceToTriangle < nodeRemoveThreshold)
				{
					nodeArray[i]->markedForDeletion = true;
					break;
				}
			}
		}
		
		//count nodes to remove
		curNode = curNode->getFirst();
		unsigned int nodesToRemove = 0;

		while(curNode->nextNode != NULL)
		{
			if(curNode->markedForDeletion)
			{
				nodesToRemove++;
			}
			curNode = curNode->nextNode;
		}

		nodeCount -= nodesToRemove;

		//remove all nodes marked for deletion
		curNode->removeMarked();

		timeStop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(timeStop - timeStart); 

		cout << "Removed " << nodesToRemove << " nodes." << endl;
		cout << "Remaining Nodes: " << nodeCount << endl;
		cout << "Node Removal Time: ";
    	outputDuration(duration.count());
    	cout << endl;
	}

	/******************************************************************
	**** .NODE FILE GENERATION
	******************************************************************/
	
	//create new empty file
	std::ofstream nodeFile;
	nodeFile.open(outputDir + "/" + fileName + "a.node", std::ofstream::out | std::ofstream::trunc);
	nodeFile.close();

	//open node file in append mode
	nodeFile.open(outputDir + "/" + fileName + "a.node", std::ios_base::app);

	//initialize node file
	nodeFile << "# Node count, 3 dim, no attribute, no boundary marker" << endl;
	nodeFile << nodeCount << " 3 0 0" << endl;
	nodeFile << "# Node index, node coordinates" << endl;

	//write node to file
	int curNodeCount = 1;
	curNode = curNode->getFirst();
	while(curNode->nextNode != NULL)
	{
		nodeFile << curNodeCount << " " << curNode->p.x << " " << curNode->p.y << " " << curNode->p.z << endl; 
		curNode = curNode->nextNode;
		curNodeCount++;
	}

	//write final node to file
	nodeFile << curNodeCount << " " << curNode->p.x << " " << curNode->p.y << " " << curNode->p.z << endl; 

	nodeFile.close();
}

/******************************************************************
* Check if a certain file exists in the current file system
* Input: 	name 	name of file to check
* Output:
* Return: 	bool 	true: file exists; false: file does not exist
*******************************************************************/
bool fileExists(const char* name)
{
	struct stat buffer;
	return (stat (name, &buffer) == 0);
}

/******************************************************************
* Print help for program
* Input:
* Output:
* Return:
*******************************************************************/
void printHelp(char* executableName)
{
	printf(HELP, executableName);
}

/*----------------------------------------------------------------*/
/* MAIN ROUTINE */
/*----------------------------------------------------------------*/

int main(int argc, char *argv[]) 
{
	//validate and read parameters
	char* objFile = argv[1];
	int gridNodeCount = 0;
	bool honeycombMode = 0;
	bool improve = 0;

	for(int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-?"))        //show help
        {
            printf(HELP, argv[0]);
            return -1;
        }
        else if (!strcmp(argv[i], "-c"))
        {
            honeycombMode = 1;
        }
        else if (!strcmp(argv[i], "-i"))
        {
            improve = 1;
        }
        else if (!strcmp(argv[i], "-f"))
        {
            i++;
            if(i < argc)
            {
                objFile = argv[i];
            }
            else
            {
                printf(HELP, argv[0]);
                return -1;
            }
        }
        else if (!strcmp(argv[i], "-n"))
        {
            i++;
            if(i < argc)
            {
                try
				{
					gridNodeCount = stoi(argv[i]);
				}
				catch(invalid_argument const &e)
				{
					printf(HELP, argv[0]);
					throw std::invalid_argument(" Error: Invalid Argument for parameter -n");
				}
				catch(out_of_range const &e)
				{
					printf(HELP, argv[0]);
					throw std::out_of_range(" Error: Integer Overflow for parameter -n");
				}
				if(gridNodeCount <= 0)
				{
					printf(HELP, argv[0]);
					throw std::out_of_range(" Error: Node Count must be >0 for parameter -n");
				}
            }
            else
            {
                printf(HELP, argv[0]);
                return -1;
            }
        }
	}

	string objFileString(objFile);
	string workingDir;
	string fileName;

	//check if file exists and extract directory & file name
	if(fileExists(objFile))
	{
		if(objFile[0] != '.' && objFile[0] != '/')
		{
			printHelp(argv[0]);
			throw std::runtime_error(" Error: filePath must be a valid file path starting with './' or '/'");
			return 0;
		}
		const size_t lastSlash = objFileString.rfind('/');
		workingDir = objFileString.substr(0, lastSlash+1);
		const size_t fileExtension = objFileString.rfind('.');
		fileName = objFileString.substr(lastSlash+1, fileExtension-1);
	}
	else
	{
		printHelp(argv[0]);
		throw std::runtime_error(" Error: File '" + objFileString + "' does not exist");
		return 0;
	}

	cout << "Info: Loading .obj file '" << objFile << "' ..." << endl;
	
	//load .obj file
 	OBJ obj;
 	if(!obj.loadObjFile(objFile))
 	{
 		cout << "Error: Loading .obj file '" << objFile << "' failed" << endl;
 		return -1;
 	}
	
 	cout << endl << "Info: Finished loading .obj file..." << endl;

 	cout << endl << "Info: Generating Grid..." << endl;

 	//generate the grid
 	generateNodeGrid(obj, gridNodeCount, workingDir, fileName, honeycombMode, improve);
 	cout << endl << "Info: Finished!" << endl;

    return 0;
}
