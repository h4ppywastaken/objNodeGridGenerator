# objNodeGridGenerator
objNodeGridGenerator is a software written in C++. Generates a node grid inside an Axis-Aligned Bounding Box for a given .obj input file and node amount. Saves the nodes to a TetGen (http://www.wias-berlin.de/software/index.jsp?id=TetGen) .node file format.

## Prerequisites

 - gcc/make (for building from source with [Makefile](makefile))

## Build

Use "make" with [Makefile](Makefile) or use gcc (see [Makefile](makefile) for gcc commands).

## Usage

Execute in console from within objNodeGridGenerator folder with command ".\objNodeGridGenerator [parameters]".

## Parameters

Use ".\objNodeGridGenerator -?" or ".\objNodeGridGenerator -h" to see the help context.

## Main Contributors

[h4ppywastaken](https://github.com/h4ppywastaken)

## Contributing

Feel free to open an issue or create a pull request at any time.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) file for details.

## Other Libraries Included

 - [Bly7/OBJ-Loader](https://github.com/Bly7/OBJ-Loader)
 
## Credit

Credit to [Bly7](https://github.com/Bly7) for developing and providing a .obj loader library under the MIT License.
