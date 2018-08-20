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

    // A proper implementation of a default constructor that initializes 
    // the matrix as an identity matrix of a given dimension.
    public Matrix(int dimension) {
        rows = dimension;
        columns = dimension;

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
        matrix = Vector.transpose(list, dimension);
    }

    public Matrix(Matrix m) {
        this.columns = m.getNumCols();
        this.rows = m.getNumRows();

        List<Vector> newVector = new ArrayList<Vector>();
        for(int i = 0; i < m.getNumRows(); i++) {
            Vector v = new Vector(m.getVectorAtIndex(i));
            newVector.add(v);
        }
        this.matrix = newVector;
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

        tempMatrix = Vector.transpose(tempMatrix, other.getNumCols());

        return new Matrix(tempMatrix, rows);
    }

    // An implementation of a function that performs Gauss-Jordan 
    // Elimination to find the determinant of the matrix.
    public Double det () {
        if(rows != columns) {
            return null;
        }

        Matrix m = new Matrix(this);
        double det = 1.0;
        int span = 0;

        det = m.rowEchelon(columns, det);
        det = m.reducedRowEchelon(columns, det);
        for(int i = 0; i < matrix.size(); i++) {
            for(int j = i; j < rows; j++){
                if(m.getMatrix().get(i).getDataAtIndex(j) != 0) {
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
        

        if(det() == null ||det() == 0.0) {
            return null;
        }
        
        Matrix m = new Matrix(rows);
        List<Vector> vectors = new ArrayList<Vector>(rows);

        for(int i = 0; i < rows; i++) {
            double[] tempRow = new double[rows * 2];

            for(int j = 0; j < rows * 2; j++) {
                if(j < rows) {
                    tempRow[j] = matrix.get(i).getDataAtIndex(j % rows);
                } else {
                    tempRow[j] = m.getVectorAtIndex(i).getDataAtIndex(j % rows);
                }
            }

            vectors.add(new Vector(tempRow, rows));
        }

        vectors = Vector.transpose(vectors, rows * 2);

        Matrix reduceInverse = new Matrix(vectors, rows);
        double det = reduceInverse.rowEchelon(reduceInverse.getNumRows(), 1);
        reduceInverse.reducedRowEchelon(reduceInverse.getNumRows(), det);

        List<Vector> vectors2 = new ArrayList<Vector>(rows);

        for(int i = 0; i < reduceInverse.getNumRows(); i++) {
            double[] tempRow = new double[reduceInverse.getNumRows()];

            for(int j = 0; j < reduceInverse.getNumCols(); j++) {
                if(j >= reduceInverse.getNumRows()) {
                    tempRow[j % rows] = reduceInverse.getVectorAtIndex(i).getDataAtIndex(j);
                } 
            }

            vectors2.add(new Vector(tempRow, rows));
        }

        vectors2 = Vector.transpose(vectors2, rows);

        return new Matrix(vectors2, rows);
    }

    public Matrix transpose() {
        return new Matrix(matrix, columns);
    }

    public int getNumRows() {
        return rows;
    }

    public int getNumCols() {
        return columns;
    }

    public List<Vector> getMatrix(){
        return matrix;
    }

    public Vector getVectorAtIndex(int i) {
        return matrix.get(i);
    }

    public double rowEchelon(int columns, double determinant){
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

				Vector base = new Vector(matrix.get(startingRow));
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

	public double reducedRowEchelon(int columns, double determinant){
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
        // Vector v = new Vector(new double[]{1,2}, dim);
        // Vector v2 = new Vector(new double[]{1,3}, dim);
        // Vector v3 = new Vector(new double[]{1,2}, dim);

        // List<Vector> list = new ArrayList<>(3);
        // list.add(v);
        // list.add(v2);
        // list.add(v3);
        // Matrix matrix = new Matrix(list, dim);
        Vector v = new Vector(new double[]{4,1}, 2);
        Vector v2 = new Vector(new double[]{6,5}, 2);
        Vector v3 = new Vector(new double[]{2,4}, 2);
  
        List<Vector> list = new ArrayList<>(dim);
        list.add(v);
        list.add(v2);
        list.add(v3);
        Matrix matrix = new Matrix(list, 2);


        // ----------- FOR MULTIPLICATION TESTING
        Vector v4 = new Vector(new double[]{1, 0, 0}, 3);
        Vector v5 = new Vector(new double[]{0, 2, 0}, 3);
        Vector v6 = new Vector(new double[]{0, 0, 2}, dim);

        List<Vector> list2 = new ArrayList<>(3);
        list2.add(v4);
        list2.add(v5);
       list2.add(v6);
        Matrix matrix2 = new Matrix(list2, 3);

        Matrix mul = matrix.times(matrix2);

        System.out.println("------ MULTIPLICATION ------");
        if(mul != null) {
            for(int i = 0; i < mul.getNumRows(); i++) {
                for(int j = 0; j < mul.getNumCols(); j++) {
                    System.out.print(mul.getVectorAtIndex(i).getDataAtIndex(j) + "\t");
                }
    
                System.out.println();
            }
        } else {
            System.out.println("null");
        }

        System.out.println();

        // ----------- FOR INVERSE TESTING
        Matrix inverse = matrix.inverse();

        System.out.println("------ INVERSE ------");
        if(inverse != null) {
            for(int i = 0; i < inverse.getNumRows(); i++) {
                for(int j = 0; j < inverse.getNumRows(); j++) {
                    System.out.print(inverse.getVectorAtIndex(i).getDataAtIndex(j) + "\t");
                }
    
                System.out.println();
            }
        } else {
            System.out.println("null");
        }

        System.out.println();

        System.out.println("------ DETERMINANT ------");
        // ----------- FOR DETERMINANT TESTING
        System.out.println(matrix.det());

        System.out.println();

        System.out.println("------ TRANSPOSE ------");
        // ----------- FOR TRANSPOSE TESTING
        Matrix transpose = matrix.transpose();
        for(int i = 0; i < transpose.getNumRows(); i++) {
            for(int j = 0; j < transpose.getNumCols(); j++) {
                System.out.print(transpose.getVectorAtIndex(i).getDataAtIndex(j) + "\t");
            }

            System.out.println();
        }
    }
}
