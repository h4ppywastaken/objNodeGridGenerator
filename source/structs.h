namespace structs
{
	/*------------------------------------------------------------------
	| Struct for standard Vector operations in 3D 
	| (used for points, vectors, and colors)
	------------------------------------------------------------------*/
	struct Vector 
	{        
	    double x, y, z;           /* Position XYZ or color RGB */

	    Vector(const Vector &b) : x(b.x), y(b.y), z(b.z) {}
	    Vector(const Vector3 &b) : x(b.X), y(b.Y), z(b.Z) {}
	    Vector(double x_=0, double y_=0, double z_=0) : x(x_), y(y_), z(z_) {}
	    
	    bool operator==(const Vector &b) const 
	    {
	        return (x == b.x) && (y == b.y) && (z == b.z);
	    }

	    Vector operator+(const Vector &b) const 
	    {
	        return Vector(x + b.x, y + b.y, z + b.z);
	    }

	    Vector operator-(const Vector &b) const
	    {
	        return Vector(x - b.x, y - b.y, z - b.z);
	    }

	    Vector operator/(double c) const 
	    {
	        return Vector(x / c, y / c, z / c);
	    }

	    Vector operator*(double c) const 
	    {
	        return Vector(x * c, y * c, z * c);
	    }

	    friend Vector operator*(double c, const Vector &b) 
	    { 
	        return b * c; 
	    }
	    
	    Vector operator%(Vector &b)
	    {
			return Vector(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
		}

        bool operator<(Vector &b)
	    {
			return  ((x+y+z) < (b.x + b.y + b.z));
		}

	    const double LengthSquared() const
	    {
	        return x * x + y * y + z * z;
	    }

	    const double Length() const
	    {
	        return sqrt(LengthSquared());
	    }

	    Vector MultComponents(const Vector &b) const
	    {
	        return Vector(x * b.x, y * b.y, z * b.z);
	    }

	    const Vector Normalized() const
	    {
	        return Vector(x, y, z) / sqrt(x*x + y*y + z*z);
	    }

	    const double Dot(const Vector &b) const 
	    {
	        return x * b.x + y * b.y + z * b.z;
	    }

	    const Vector Cross(const Vector &b) const
	    {
	        return Vector((y * b.z) - (z * b.y), 
	                      (z * b.x) - (x * b.z), 
	                      (x * b.y) - (y * b.x));
	    }
	     
	    const double Max() 
	    {
	        return fmax(z, fmax(x, y));
	    }

        const double Min() 
	    {
	        return fmin(z, fmin(x, y));
	    }

	    Vector& clamp() 
	    {
	        x = x<0 ? 0.0 : x>1.0 ? 1.0 : x;
	        y = y<0 ? 0.0 : y>1.0 ? 1.0 : y;
	        z = z<0 ? 0.0 : z>1.0 ? 1.0 : z;
	        return *this;   
	    }

	    Vector& clampScale()
	    {
	    	double maxRGB = fmax(z, fmax(x, y));;
	    	if(maxRGB > 1.0)
	    	{
	    		x = x / maxRGB;
	    		y = y / maxRGB;
	    		z = z / maxRGB;
	    	}
	        return *this;  
	    }

	    const void printVector() const
	    {
	        cout << "\t x " << x << "\t y " << y << "\t z " << z << endl;
	    }

	    //Geo Servant Format: https://www.matheretter.de/geoservant
	    const void printGeoServantVector() const
	    {
	        cout << "vector(0|0|0 " << x << "|" << y << "|" << z << ")" << endl;
	    }

	    //Geo Servant Format: https://www.matheretter.de/geoservant
	    const void printGeoServantPoint() const
	    {
	        cout << "point(" << x << "|" << y << "|" << z << ")" << endl;
	    }
	};


	/*------------------------------------------------------------------
	| Struct for edges with start and end vector 
	| (used for edges in triangles)
	| Based on https://stackoverflow.com/questions/2924795/fastest-way-to-compute-point-to-triangle-distance-in-3d
	------------------------------------------------------------------*/
	struct Edge
	{
	    Vector a;
	    Vector b;
	    Vector delta;

	    Edge(const Vector &a_, const Vector &b_) : a(a_), b(b_)
	    {
	        delta = b - a;
	    }

	    Edge(){}

	    const Vector PointAt(double t) const
	    {
			return a + t * delta;
	    }

	    const double LengthSquared() const
	    {
	    	return delta.LengthSquared();
	    }

	    //projects vector onto the edge
	    const double Project(Vector p) const
	    {
	    	return (p - a).Dot(delta) / LengthSquared();
	    }
	};

	/*------------------------------------------------------------------
	| Struct for planes defined by a point and a vector direction
	| (used for edge planes in triangles)
	| Based on https://stackoverflow.com/questions/2924795/fastest-way-to-compute-point-to-triangle-distance-in-3d
	------------------------------------------------------------------*/
	struct Plane
	{
	    Vector point;
	    Vector normal;

	    Plane(const Vector &point_, const Vector &normal_) : point(point_)
	    {
	    	normal = normal_.Normalized();
	    }
	    Plane(){}

	    //checks if plane is above a vector
	    const bool isAbove(Vector p) const
	    {
	    	return normal.Dot(p - point) < 0;
	    }

	    //projects vector onto the plane
	    const Vector Project(Vector p) const
	    {
	    	Vector v = p - point;
	    	double d = v.Dot(normal);
	    	return p - d*normal;
	    }
	};

	/*------------------------------------------------------------------
	| Helper function to calculate mid-point of 3 given vectors
	------------------------------------------------------------------*/
	const Vector getMiddle(Vector p0, Vector p1, Vector p2)
	{
	    return (p0 * 1 / 3) + (p1 * 1 / 3) + (p2 * 1 / 3);
	}

	/*------------------------------------------------------------------
	| Helper function to compare two double values with a certain error tolerance
	------------------------------------------------------------------*/
	const bool withinDoubleError(double value1, double value2)
	{
		double doubleError = 0.00000001;
		return value1 < (value2 + doubleError) && value1 > (value2 - doubleError);
	}

	/*------------------------------------------------------------------
	| Basic geometric element of scene description;
	------------------------------------------------------------------*/
	struct Triangle
	{
	    Vector p0, p1, p2, middle;
	    Vector normal;
	    Edge edge01, edge12, edge20;
	    Plane plane01, plane12, plane20;

	    Triangle(const Vector &p0_, const Vector &p1_, const Vector &p2_) : p0(p0_), p1(p1_), p2(p2_)
	    {
	    	edge01 = Edge(p0, p1);
	        edge12 = Edge(p1, p2);
	        edge20 = Edge(p2, p0);

	        normal = ((edge01.delta).Cross(edge20.delta)).Normalized();
	        middle = getMiddle(p0,p1,p2);
	        
	        plane01 = Plane(p0, normal.Cross(edge01.delta));
    		plane12 = Plane(p1, normal.Cross(edge12.delta));
    		plane20 = Plane(p2, normal.Cross(edge20.delta));

    		//Sanity Checks for flipped normals
    		if(!plane01.isAbove(p2)) plane01.normal = (-1) * plane01.normal;
    		if(!plane12.isAbove(p0)) plane12.normal = (-1) * plane12.normal;
    		if(!plane20.isAbove(p1)) plane20.normal = (-1) * plane20.normal;
	    }
	    
	    Triangle(){}

	    //Geo Servant Format: https://www.matheretter.de/geoservant
	    const void printGeoServantTriangle() const
	    {
	        cout << "triangle(" << p0.x << "|" << p0.y << "|" << p0.z << " " << p1.x << "|" << p1.y << "|" << p1.z << " " << p2.x << "|" << p2.y << "|" << p2.z << ")" << endl;
	    }

	    //Checks if a vector is pointing inside the triangle
	    const bool CheckCollision(Vector p)
	    {
	    	//point is corner vertex of triangle
	    	if(p == p0 || p == p1 || p == p2)
				return 1;

			double planeColl = normal.Dot(p - middle);

			//point not within triangle plane
			if(!withinDoubleError(planeColl, 0))
			{
				return 0;
			}

			double pedge01 = edge01.Project(p);
			double pedge12 = edge12.Project(p);
	        double pedge20 = edge20.Project(p);
	        
	        //point is on triangle edge
	        if((withinDoubleError(plane01.normal.Dot(p-p0), 0) && 0 < pedge01 && pedge01 < 1)
				|| (withinDoubleError(plane12.normal.Dot(p-p1), 0) && 0 < pedge12 && pedge12 < 1)
				|| (withinDoubleError(plane20.normal.Dot(p-p2), 0) && 0 < pedge20 && pedge20 < 1))
			{
				return 1;
			}

	        //point within triangle
			if(plane01.isAbove(p) && plane12.isAbove(p) && plane20.isAbove(p))
			{
				return 1;
			}

			return 0;
	    }

	    //Gets Closest point in the triangle in respect to a vector
	    //Based on https://stackoverflow.com/questions/2924795/fastest-way-to-compute-point-to-triangle-distance-in-3d
	    const Vector ClosestPoint(Vector p)
	    {
	    	//check if point is closest to a triangle corner vertex
	        double pedge01 = edge01.Project(p);
	        double pedge20 = edge20.Project(p);

	        if (pedge20 >= 1 && pedge01 <= 0)
	            return p0;

	        double pedge12 = edge12.Project(p);

	        if (pedge01 >= 1 && pedge12 <= 0)
	            return p1;

	        if (pedge12 >= 1 && pedge20 <= 0)
	            return p2;

	        //check if point is next to the triangle "above" an edge
	        if ((0.0 < pedge01 && pedge01 < 1.0) && !plane01.isAbove(p))
	            return edge01.PointAt(pedge01);

	        if ((0.0 < pedge12 && pedge12 < 1.0) && !plane12.isAbove(p))
	            return edge12.PointAt(pedge12);

	        if ((0.0 < pedge20 && pedge20 < 1.0) && !plane20.isAbove(p))
	            return edge20.PointAt(pedge20);

	        //point is above or below triangle
	        return Plane(middle, normal).Project(p);
	    }
	};


	/*------------------------------------------------------------------
	| Struct for Nodes
	| (double-linked list)
	------------------------------------------------------------------*/
	struct Node
	{
		Vector p;
		Node *nextNode, *previousNode;
		bool markedForDeletion;

		Node(const Vector &p_) : p(p_)
	    {
	    	previousNode = NULL;
	    	nextNode = NULL;
	    	markedForDeletion = false;
	    }

	    Node()
	    {
	    	previousNode = NULL;
	    	nextNode = NULL;
	    	markedForDeletion = true;
	    }

	    Node * appendNode(Node *node)
		{
			if(markedForDeletion)
			{
				(*this) = Node(node->p);
				return this;
			}
			else
			{
				nextNode = node;
				node->previousNode = this;
				return node;
			}
		}

		Node * prependNode(Node *node)
		{
			if(markedForDeletion)
			{
				(*this) = Node(node->p);
				return this;
			}
			else
			{
				previousNode = node;
				node->nextNode = this;
				return node;
			}
		}

		Node * getFirst()
		{
			Node *curNode = this;
			while(curNode->previousNode != NULL)
			{
				curNode = curNode->previousNode;
			}
			return curNode;
		}

		Node * removeNode()
		{
			Node *temp = this;
			Node *returnNode;
			if(previousNode != NULL && nextNode != NULL)
			{
				previousNode->nextNode = nextNode;
				nextNode->previousNode = previousNode;
				returnNode = nextNode;
			}
			else if(nextNode != NULL)
			{
				nextNode->previousNode = NULL;
				returnNode = nextNode;
			}
			else if(previousNode != NULL)
			{
				previousNode->nextNode = NULL;
				returnNode = previousNode;
			}
			else
			{
				returnNode = new Node();
			}
			//free(temp);
			return returnNode;
		}

		//Removes all nodes from the list that have markedForDeletion == true
		Node * removeMarked()
		{
			Node *curNode = getFirst();
			if(curNode->markedForDeletion)
			{
				curNode->removeNode()->removeMarked();
			}
			while(curNode->nextNode != NULL)
			{
				if(curNode->markedForDeletion)
				{
					curNode = curNode->removeNode();
				}
				else
				{
					curNode = curNode->nextNode;
				}
			}
			return getFirst();
		}
	};

	/*------------------------------------------------------------------
	| Struct for .obj file loading with Circumscribing Cuboid calculation
	------------------------------------------------------------------*/
	struct OBJ
	{
		//all object triangles
	    Triangle *triangles;
	    unsigned int triangleCount;

	    //Circumscribed rectangular cuboid of the object
	    double minX, minY, minZ, maxX, maxY, maxZ;
	    double edgeX, edgeY, edgeZ;

	    OBJ(){};

	    bool loadObjFile(char *objfile)
		{
			// Initialize Loader
			objl::Loader Loader;

		    bool loadout = Loader.LoadFile(objfile);

			objl::Mesh mesh = Loader.LoadedMeshes[0];

			triangleCount = mesh.Indices.size()/3;

		    triangles = new Triangle[triangleCount];

		    if(!loadout)
		    {
		    	cout << "Error loading obj file:" << objfile << endl << flush;
		    	return 0;
		    }

			cout << "Loaded Mesh: " << mesh.MeshName << endl << flush;

			//initialize circumscribed rectangular cuboid
			minX = mesh.Vertices[mesh.Indices[0]].Position.X;
			minY = mesh.Vertices[mesh.Indices[0]].Position.Y;
			minZ = mesh.Vertices[mesh.Indices[0]].Position.Z;

			maxX = minX;
			maxY = minY;
			maxZ = minZ;

			//Save all triangles in the mesh and calculate normal & middle of each triangle (see Triangle struct)
			for (unsigned int j = 0; j < triangleCount*3; j += 3)
			{
		        Vector3 p0 = mesh.Vertices[mesh.Indices[j]].Position;
		        Vector3 p1 = mesh.Vertices[mesh.Indices[j + 1]].Position;
		        Vector3 p2 = mesh.Vertices[mesh.Indices[j + 2]].Position;
		        
				triangles[j/3] = 
				Triangle(
					Vector(p0),
					Vector(p1),
					Vector(p2)
					);

				//determine circumscribed rectangular cuboid
				if(p0.X < minX) minX = p0.X;
				if(p0.Y < minY) minY = p0.Y;
				if(p0.Z < minZ) minZ = p0.Z;

				if(p0.X > maxX) maxX = p0.X;
				if(p0.Y > maxY) maxY = p0.Y;
				if(p0.Z > maxZ) maxZ = p0.Z;

				if(p1.X < minX) minX = p1.X;
				if(p1.Y < minY) minY = p1.Y;
				if(p1.Z < minZ) minZ = p1.Z;

				if(p1.X > maxX) maxX = p1.X;
				if(p1.Y > maxY) maxY = p1.Y;
				if(p1.Z > maxZ) maxZ = p1.Z;

				if(p2.X < minX) minX = p2.X;
				if(p2.Y < minY) minY = p2.Y;
				if(p2.Z < minZ) minZ = p2.Z;

				if(p2.X > maxX) maxX = p2.X;
				if(p2.Y > maxY) maxY = p2.Y;
				if(p2.Z > maxZ) maxZ = p2.Z;
			}

			//calculate edge length of circumscribed rectangular cuboid
			edgeX = fabs(maxX - minX);
			edgeY = fabs(maxY - minY);
			edgeZ = fabs(maxZ - minZ); 

			return 1;
		}
	};
}