/**
 * @author BATOSALEM, Angelika && CORTEZ, Louise && RIVERA, Sophia
 * @section S17
 */

import java.util.ArrayList;
import java.util.List;

public class Matrix {
    // The usage of an Array/List-like structure to store Matrix data 
    // as a list of Vectors. You may also store the Matrix as a 2d array.
    private List<Vector> matrix;

    // The usage of immutable Integer variables to hold values for 
    // the number of rows/columns.
    private final int rows;
    private final int columns;
    
    private int dimension; // ?? not sure

    // A proper implementation of a default constructor that initializes 
    // the matrix as an identity matrix of a given dimension.
    public Matrix(int dimension) {
        rows = dimension;
        columns = dimension;
        this.dimension = dimension;
        double[] row = new double[dimension];
        matrix = new ArrayList<Vector>();

        for(int i = 0; i < dimension; i++) {
            for(int j = 0; j < dimension; j++) {
                if(i == j) {
                    row[j] = 1;
                } else {
                    row[j] = 0;
                }
            }

            Vector v = new Vector(row, dimension);
            matrix.add(v);
        }
    }

    // A proper implementation of a constructor, converting an 
    // already-existing array/list of data from a rudimentary data 
    // structure into the vector class.
    public Matrix (List<Vector> list, int dimension) {
        columns = list.size();
        rows = dimension;
        this.dimension = dimension;

        // TODO remove transpose in Vector class (?)
        matrix = Vector.transpose(list, dimension);
    }

    // An implementation of function for matrix multiplication.
    public Matrix times (Matrix other) {
        // Errors for size mismatches when multiplying matrices must also be handled.
        if(columns != other.getNumRows()) {
            return null;
        }

        List<Vector> tempMatrix = new ArrayList<Vector>();
        double[][] multiplied = new double[rows][other.getNumCols()];

        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < other.getNumCols(); j++) {
                for(int k = 0; k < other.getNumRows(); k++) {
                    multiplied[i][j] += matrix.get(i).getDataAtIndex(k) * other.getVectorAtIndex(k).getDataAtIndex(j);
                }
            }

            Vector v = new Vector(multiplied[i], other.getNumCols());
            tempMatrix.add(v);
        }

        return new Matrix(tempMatrix, rows);
    }

    // An implementation of a function that performs Gauss-Jordan 
    // Elimination to find the determinant of the matrix.
    public double det () {
        if(rows != columns) {
            return 0.0;
        }

        double det = 1.0;
        int span = 0;
        det = rowEchelon(matrix, columns, det);
        det = reducedRowEchelon(matrix, columns, det);
        for(int i = 0; i < matrix.size(); i++) {
            for(int j = i; j < rows; j++){
                if(matrix.get(i).getDataAtIndex(j) != 0) {
                    span++;
                    break;
                }
            }
        }


        return (rows == span)?1/det:0;
    }

        // An implementation of a function that finds the inverse of the matrix.
        public Matrix inverse() {
            // // The function must return a null value if the matrix has no inverse.
            Matrix temp = new Matrix(matrix, rows);

            if(temp.det() == 0.0) {
                return null;
            }
            
            Matrix m = new Matrix(rows);
            List<Vector> vectors = new ArrayList<Vector>(rows);
    
            for(int i = 0; i < rows; i++) {
                double[] tempRow = new double[rows * 2];
    
                for(int j = 0; j < rows * 2; j++) {
                    if(j < rows) {
                        tempRow[j] = matrix.get(j % rows).getDataAtIndex(i);
                    } else {
                        tempRow[j] = m.getVectorAtIndex(i).getDataAtIndex(j % rows);
                    }
                }
    
                vectors.add(new Vector(tempRow, rows));
            }

            for(int i = 0; i < rows; i++) {
                for(int j = 0; j < rows * 2; j++) {
                    System.out.print(vectors.get(i).getDataAtIndex(j) + "\t");
                }

                System.out.println();
            }
            return new Matrix(vectors, rows);
        }

    public int getNumRows() {
        return rows;
    }

    public int getNumCols() {
        return columns;
    }

    public Vector getVectorAtIndex(int i) {
        return matrix.get(i);
    }

    public double rowEchelon(List<Vector> matrix, int columns, double determinant){
		int startingRow = 0;
		for (int col = 0; col < columns; col++) {
			int max = startingRow;
			
			while(matrix.get(max).getDataAtIndex(col) == 0 && max < matrix.size()-1)
				max++;
			
			if(matrix.get(max).getDataAtIndex(col) != 0) {
				if(max != startingRow) {
					Vector v = matrix.get(max);
					matrix.set(max, matrix.get(startingRow));
					matrix.set(startingRow, v);

                    determinant = -1 * determinant;
				}

				Vector base = new Vector(matrix.get(max));
				for(int row = startingRow+1; row < matrix.size(); row++) {
					if(matrix.get(row).getDataAtIndex(col) != 0) {
						double factor = -1 * matrix.get(row).getDataAtIndex(col);

                        determinant = base.getDataAtIndex(col) * determinant;
						matrix.get(row).scale(base.getDataAtIndex(col)).add(base.scale(factor));
					
						base.scale(1/factor);
					}

				}

				startingRow++;
			}
		}

        return determinant;
	}

	public double reducedRowEchelon(List<Vector> matrix, int columns, double determinant){
		for(int i = matrix.size() - 1; i >= 0 ; i--) {
			int nonzeroIndex = -1;
			for(int j = 0; j < columns; j++) {
				if(matrix.get(i).getDataAtIndex(j) != 0) {
					nonzeroIndex = j;
					break;
				}
			}

			if(nonzeroIndex != -1) {
				double factor = matrix.get(i).getDataAtIndex(nonzeroIndex);
				matrix.get(i).scale(1/factor);
                determinant = determinant * 1/factor;

				if(i > 0) {
					for(int k = i - 1; k >= 0; k--) {
						if(matrix.get(k).getDataAtIndex(nonzeroIndex) != 0) {
							Vector base = new Vector(matrix.get(i));
							double factor2 = -1 * matrix.get(k).getDataAtIndex(nonzeroIndex);

							matrix.get(k).scale(base.getDataAtIndex(nonzeroIndex)).add(base.scale(factor2));
							base.scale(1/factor2);
                            determinant = determinant * base.getDataAtIndex(nonzeroIndex);
						}
					}
				}
			}
		}

        return determinant;

	}

    // For testing
    public static void main(String[] args) {
        Matrix m = new Matrix(3);
        int dim = 3;
        /*Vector v = new Vector(new double[]{1,2}, dim);
        Vector v2 = new Vector(new double[]{1,2}, dim);*/
        Vector v = new Vector(new double[]{4,1,2}, dim);
        Vector v2 = new Vector(new double[]{6,5,7}, dim);
        Vector v3 = new Vector(new double[]{2,4,3}, dim);
  
        List<Vector> list = new ArrayList<>(dim);
        list.add(v);
        list.add(v2);
        list.add(v3);
        Matrix matrix = new Matrix(list, dim);


        // ----------- FOR MULTIPLICATION TESTING
        // Vector v4 = new Vector(new double[]{2, 0, 0}, dim);
        // Vector v5 = new Vector(new double[]{0, 2, 0}, dim);
        // Vector v6 = new Vector(new double[]{0, 0, 2}, dim);

        // List<Vector> list2 = new ArrayList<>(dim);
        // list2.add(v4);
        // list2.add(v5);
        // list2.add(v6);
        // Matrix matrix2 = new Matrix(list2, dim);

        // Matrix mul = matrix.times(matrix2);

        // for(int i = 0; i < mul.getNumRows(); i++) {
        //     for(int j = 0; j < mul.getNumCols(); j++) {
        //         System.out.print(mul.getVectorAtIndex(i).getDataAtIndex(j) + "\t");
        //     }

        //     System.out.println();
        // }

        matrix.inverse();
        // System.out.println(matrix.det());
    }
}