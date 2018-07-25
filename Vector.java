import java.util.ArrayList;
import java.util.List;
import java.util.Collections;

public class Vector {
	static int min(int a, int b) { 
        return a < b ? a : b; 
         
    }

	//A proper implementation of a vector function via the usage of a List-like data structure. (5 points)


	//  The usage of an Array/List-like structure to store the Vector data.
	private ArrayList<Double> data;
	//  The usage of an immutable Integer variable to hold a value for Vector dimension.
	private final int dimension;

	//  A proper implementation of a default constructor that initializes the vector as a zero vector of a given dimension.
	//      Constructor definition to be used: Vector (int dimension)

	public Vector(int dimension) {
		data = new ArrayList<>(Collections.nCopies(dimension, 0.0));
		this.dimension = dimension;
	}

	//        A proper implementation of a constructor, converting an already-existing array/list of data from a rudimentary data structure into the vector class.
	//            Constructor definition to be used: Vector (double[ ] array, int dimension)

	public Vector(double[] array, int dimension) {
		this.data = new ArrayList<>(dimension); // not sure if need pa?
		this.dimension = dimension;

		for(int i = 0; i < array.length; i++) {
			this.data.add(array[i]);
		}
	}

	//copying
	public Vector(Vector newVector) {
		this.data = new ArrayList<>(newVector.getVector());
		this.dimension =  data.size();
	}

	//    An implementation of functions for vector scaling and vector addition. (10 points)
	//        A proper implementation of a function for vector scaling.
	//            Function header to be used: Vector scale (int scalar)
	//  Usage example: Assuming a Vector v and int b exists, v.scale(b) should scale the elements of v by b and return the scaled vector v. The elements inside v must be changed and be correctly scaled by b.

	public Vector scale(int scalar) {

		for(int i = 0; i < data.size(); i++) {
			data.set(i,data.get(i) * scalar);
		}
		return this;
	}

	public Vector reduce(int scalar) {

		for(int i = 0; i < data.size(); i++) {
			data.set(i,data.get(i) / scalar);
		}
		return this;
	}
	//        A proper implementation of a function for vector addition. Errors for size mismatches when adding vectors must also be handled.
	//            Function header to be used: Vector add (Vector addend)
	//            Usage example: Assuming both Vector v and w exist, v.add(w) should return the vector sum between v and w. The elements inside v must be changed and be a correct result of the operation of Vector addition between v and w.

	public Vector add(Vector addend) {
		if(data.size() != addend.getSize()) {
			// error handling
			return null;
		} 


		for(int i = 0; i < data.size(); i++) {
			data.set(i,data.get(i) + addend.getDataAtIndex(i));
		}

		return this;
	}
	//    An implementation of a function that performs Gauss-Jordan Elimination on a given set of vectors. (30 points)
	
	//        The function must be static-like in nature, and must be callable from the Vector class. See usage example for more details.
	//        Function header to be used: Vector Gauss_Jordan (List<Vector> vectors, int dimension, Vector constants)
	//        The function must be a proper implementation of Gauss-Jordan Elimination, which reduces the given list of vectors into unit vectors via row operations.
	//        Usage example: Given a list of vectors vecList, an integer dim, and a Vector c, Vector.Gauss_Jordan (vecList, dim, c) should return a Vector containing the solution to the corresponding system of linear equations. Ex. [x y z w] = [2 1 3 5]
	
	public static Vector Gauss_Jordan(List<Vector> vectors, int dimension, Vector constants) {
		Vector.swap(vectors, dimension, constants);
		Vector.rowEchelon(vectors, dimension, constants);
		// Vector.reducedRowEchelon(vectors,dimension, constants);
		return constants;
	}

	public static void swap(List<Vector> vectors, int dimension, Vector constants){
		int currInd = 0, currZero = 0;
		boolean switched = true;
		while(currInd < vectors.size() && switched && currZero < dimension) {

			switched = false;
			for(int i = currInd; i < vectors.size(); i++) {
				int j = currInd;
				if(vectors.get(i).getDataAtIndex(currZero) != 0 
				&& vectors.get(i).getDataAtIndex(currZero) > vectors.get(j).getDataAtIndex(currZero)) {
					
					if(i != j) {
						Vector v = vectors.get(j);
						vectors.set(j, vectors.get(i));
						vectors.set(i, v);
						// add switch for constants
						double temp = constants.getDataAtIndex(j);
						constants.setValue(j, constants.getDataAtIndex(i));
						constants.setValue(i, temp);
					}

					j++;
					currInd = j;
					switched = true;
				}
			}
			currZero++;
		}
		//return vectors;
	}

	public static void rowEchelon(List<Vector> vectors, int dimension, Vector constants){
		// for(int i = 0; i < vectors.size() && i < dimension; i++) {
		// 	for(int j = i+1; j < vectors.size() && j < dimension; j++) {
		// 		if(vectors.get(j).getDataAtIndex(i) != 0){

		// 			Vector base = new Vector(vectors.get(i));
		// 			int factor = -1 * vectors.get(j).getDataAtIndex(i).intValue();

		// 			double first = constants.getDataAtIndex(j)*base.getDataAtIndex(i),
		// 			second = constants.getDataAtIndex(i)*factor;

		// 			vectors.get(j).scale(base.getDataAtIndex(i).intValue()).add(base.scale(factor));
		// 			//add constant part
		// 			constants.setValue(j,first+second);

		// 		}
		// 	}

		// 	System.out.println("i is: " + i);
		// 	for(int j = 0; j < vectors.size(); j ++) {
		// 		for(int k = 0; k < vectors.get(j).getSize(); k++) {
		// 			System.out.printf("%.5f ",vectors.get(j).getDataAtIndex(k));
		// 		}
		// 		System.out.println();
		// 	}

		// 	System.out.println();
		// }
		// int reduceFactor = vectors.get(dimension-1).getDataAtIndex(dimension-1).intValue();
		// if(reduceFactor != 0){
		// 	vectors.get(dimension-1).reduce(reduceFactor);
		// 	// add constant part
		// 	constants.setValue(dimension-1,constants.getDataAtIndex(dimension-1)/reduceFactor);
		// }
		// // for(int i = 0; i < vectors.size(); i++)
		// // 	if(vectors.get(i).getDataAtIndex(i).intValue() != 0)
		// // 		vectors.get(i).reduce(vectors.get(i).getDataAtIndex(i).intValue());
		// // return vectors;

		// PIVOT ALGORITHM FROM WEBSITE
		// 1. Assuming the list is sorted, start with the pivot index of the first row (first nonzero digit)
		// 2. Convert all digits in the column to zero besides the pivot index
		// 3. Convert digit in pivot index to 1 by dividing it by itself
		// 4. Move to the next row, then choose the next pivot index (first nonzero digit >= current row)
		// 5. Repeat process
		// NOT SURE: when magrereduce and stuff :'(
		int pivotIndex = 0;
		for(int i = 0; i < vectors.size(); i++) {
			for(int j = 0; j < dimension; j++) {
				if(vectors.get(i).getDataAtIndex(j) != 0) {
					pivotIndex = j;
					System.out.println("Hello the pivot is at index: " + i + " " + j);
					break;
				}
			}

			for(int k = 0; k < vectors.size(); k++) {
				if(i != k) {
					if(vectors.get(k).getDataAtIndex(pivotIndex) != 0){
						Vector base = new Vector(vectors.get(i));
						int factor = -1 * vectors.get(k).getDataAtIndex(pivotIndex).intValue();
	
						double first = constants.getDataAtIndex(k)*base.getDataAtIndex(pivotIndex),
						second = constants.getDataAtIndex(pivotIndex)*factor;
	
						vectors.get(k).scale(base.getDataAtIndex(l).intValue()).add(base.scale(factor));
						//add constant part
						constants.setValue(k,first+second);
					}
				}
			}

			int reduceFactor = vectors.get(i).getDataAtIndex(pivotIndex).intValue();
			if(reduceFactor != 0){
				vectors.get(i).reduce(reduceFactor);
				// add constant part
				constants.setValue(i,constants.getDataAtIndex(pivotIndex)/reduceFactor);
			}
		}
	}

	public static void reducedRowEchelon(List<Vector> vectors, int dimension, Vector constants){
		for(int i = min(vectors.size()-1,dimension-1); i > 0; i--) {
			
			if(vectors.get(i).getDataAtIndex(i) != 0) {
				for(int j = i-1; j >= 0; j--) {
					if(vectors.get(j).getDataAtIndex(i) != 0){
						Vector base = new Vector(vectors.get(i));
						int factor = -1 * vectors.get(j).getDataAtIndex(i).intValue();

						double first = constants.getDataAtIndex(j)*base.getDataAtIndex(i),
						second = constants.getDataAtIndex(i)*factor;

						vectors.get(j).scale(base.getDataAtIndex(i).intValue()).add(base.scale(factor));
						// add constant part
						
						constants.setValue(j,first+second);
					}
				}
			}
		}
		for(int i = 0; i < vectors.size(); i++) {
			for(int j = 0; j < dimension; j++) {
				int reduceFactor = vectors.get(i).getDataAtIndex(j).intValue();
				if(reduceFactor != 0){
					vectors.get(i).reduce(reduceFactor);
					// add constant part
					constants.setValue(i,constants.getDataAtIndex(i)/reduceFactor);
				}
			}
		}
		//return vectors;
	}
	
	//    An implementation of a function that calculates the span of a list of vectors. (5 points)
	//        The function must be static-like in nature, and must be callable from the Vector class. See usage example for more details.
	//        Function header to be used: int span (List<Vector> vectors, int dimension)
	
	public static int span(List<Vector> vectors, int dimension) {
		//  The function must call the Gauss_Jordan function within it; i.e. the calculation for span must include Gauss-Jordan Elimination.
		//        Usage example: Given a list of Vectors vecList, and an integer dim denoting the dimension of a vector inside vecList, Vector.span(vecList, dim) should return an integer variable containing the span of the set of vectors.
		int span = 0;

		Vector zeroConstants = new Vector(vectors.size());
		Vector gaussJordan = Gauss_Jordan(vectors, dimension, zeroConstants);

		for(int i = 0; i < vectors.size(); i++) {
			for(int j = i; j < dimension; j++){
				if(vectors.get(i).getDataAtIndex(j) != 0) {
					span++;
					break;
				}
			}
		}
		
		return span;
	}

	public List<Vector> transpose(List<Vector> vectors, int dimension) {
		double[][] transposedDigits = new double[dimension][vectors.size()];
		List<Vector> transposedMatrix = new ArrayList<Vector>(vectors.size());

		for(int i = 0; i < dimension; i++) {
			for(int j = 0; j < vectors.size(); j++) {
				transposedDigits[i][j] = vectors.get(j).getDataAtIndex(i);
			}
		}

		for(int i = 0; i < dimension; i++) {
			Vector v = new Vector(transposedDigits[i], vectors.size());
			transposedMatrix.add(v);
		}

		return transposedMatrix;
	}

	// for testing
	public ArrayList<Double> getVector() {
		return data;
	}

	public int getSize() {
		return data.size();
	}

	public void setValue(int index, double newValue) {
		data.set(index, newValue);
	}

	// for testing
	public Double getDataAtIndex(int i) {
		return data.get(i);
	}

	// for testing (to be removed)
	public static void main(String[] args) {
		// double[] vector = new double[]{0, 0, 0, 6, 0};
		// double[] vector2 = new double[]{0, 0, 10, 13, 1};
		// double[] vector4 = new double[]{4, 4.3, 10, 13, 9};
		// double[] vector5 = new double[]{4, 0, 10, 13, 7};
		// double[] vector3 = new double[]{0, 5, 10, 13, 3};
		int dimension = 2;
		// double[] vector = new double[]{2, 4, -2, 8, 4};
		// double[] vector2 = new double[]{3, 6, 1, 12, -2};
		// double[] vector3 = new double[]{9, 18, 1, 36, 38};
		// double[] vector4 = new double[]{2, -1, 2,-1};
		// double[] cons = new double[]{1, -2, 4};

		double[] vector = new double[]{1, 2};
		double[] vector2 = new double[]{1, 1};
		double[] vector3 = new double[]{2, 3};
		double[] cons = new double[]{3, 1, 1};

		// double[] vector = new double[]{0.02, 0.01, 0, 0};
		// double[] vector2 = new double[]{1, 2, 1, 0};
		// double[] vector3 = new double[]{0, 1, 2, 1};
		// double[] vector4 = new double[]{0, 0, 100, 200};
		// double[] cons = new double[]{0.02, 1, 4, 800};

		Vector v = new Vector(vector, dimension);
		Vector v2 = new Vector(vector2, dimension);
		Vector v3 = new Vector(vector3, dimension);
		// Vector v4 = new Vector(vector4, dimension);
		Vector constants = new Vector(cons, dimension);

		List<Vector> list = new ArrayList<>(dimension);
		list.add(v);
		list.add(v2);
		list.add(v3);
		// list.add(v4);

		/*double[] vector = new double[]{2,3};
		double[] vector2 = new double[]{1,5};
		double[] cons = new double[]{7,7};
		Vector v = new Vector(vector, dimension);
		Vector v2 = new Vector(vector2, dimension);
		Vector constants = new Vector(cons, dimension);

		List<Vector> list = new ArrayList<>(dimension);
		list.add(v);
		list.add(v2);*/


		Vector result = Vector.Gauss_Jordan(list, dimension, constants);
		Vector.span(list, dimension);
		for(int j = 0; j < list.size(); j ++) {
			for(int i = 0; i < list.get(j).getSize(); i++) {
				System.out.printf("%.5f ",list.get(j).getDataAtIndex(i));
			}
			System.out.println();
		}
		for(int i = 0; i < result.getSize(); i++) {
			System.out.print(result.getDataAtIndex(i) + " ");
		}
	}
}
