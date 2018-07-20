import java.util.ArrayList;
import java.util.List;

public class Vector {


	//A proper implementation of a vector function via the usage of a List-like data structure. (5 points)


	//  The usage of an Array/List-like structure to store the Vector data.
	private ArrayList<Integer> data;
	//  The usage of an immutable Integer variable to hold a value for Vector dimension.
	private int dimension;

	//  A proper implementation of a default constructor that initializes the vector as a zero vector of a given dimension.
	//      Constructor definition to be used: Vector (int dimension)

	public Vector(int dimension) {
		data = new ArrayList<>();
		this.dimension = dimension;
	}

	//        A proper implementation of a constructor, converting an already-existing array/list of data from a rudimentary data structure into the vector class.
	//            Constructor definition to be used: Vector (double[ ] array, int dimension)

	public Vector(double[] array, int dimension) {
		
	}

	//    An implementation of functions for vector scaling and vector addition. (10 points)
	//        A proper implementation of a function for vector scaling.
	//            Function header to be used: Vector scale (int scalar)
	//  Usage example: Assuming a Vector v and int b exists, v.scale(b) should scale the elements of v by b and return the scaled vector v. The elements inside v must be changed and be correctly scaled by b.

	public Vector scale(int scalar) {
		Vector v = new Vector(1);
		return v;
	}
	//        A proper implementation of a function for vector addition. Errors for size mismatches when adding vectors must also be handled.
	//            Function header to be used: Vector add (Vector addend)
	//            Usage example: Assuming both Vector v and w exist, v.add(w) should return the vector sum between v and w. The elements inside v must be changed and be a correct result of the operation of Vector addition between v and w.

	public Vector add(Vector addend) {
		Vector v = new Vector(1);
		return v;
	}
	//    An implementation of a function that performs Gauss-Jordan Elimination on a given set of vectors. (30 points)
	
	//        The function must be static-like in nature, and must be callable from the Vector class. See usage example for more details.
	//        Function header to be used: Vector Gauss_Jordan (List<Vector> vectors, int dimension, Vector constants)
	//        The function must be a proper implementation of Gauss-Jordan Elimination, which reduces the given list of vectors into unit vectors via row operations.
	//        Usage example: Given a list of vectors vecList, an integer dim, and a Vector c, Vector.Gauss_Jordan (vecList, dim, c) should return a Vector containing the solution to the corresponding system of linear equations. Ex. [x y z w] = [2 1 3 5]
	
	public Vector Gauss_Jordan(List<Vector> vectors, int dimension, Vector constants) {
	
		Vector v = new Vector(1);
		return v;
	}
	
	//    An implementation of a function that calculates the span of a list of vectors. (5 points)
	//        The function must be static-like in nature, and must be callable from the Vector class. See usage example for more details.
	//        Function header to be used: int span (List<Vector> vectors, int dimension)
	
	public int span(List<Vector> vectors, int dimension) {
		//  The function must call the Gauss_Jordan function within it; i.e. the calculation for span must include Gauss-Jordan Elimination.
		//        Usage example: Given a list of Vectors vecList, and an integer dim denoting the dimension of a vector inside vecList, Vector.span(vecList, dim) should return an integer variable containing the span of the set of vectors.

		return 0;
	}
	
		 


}
