/**
 * @author BATOSALEM, Angelika && CORTEZ, Louise && RIVERA, Sophia
 * @section S17
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Collections;

public class Vector {
	private ArrayList<Double> data;
	private final int dimension;

	public Vector(int dimension) {
		data = new ArrayList<>(Collections.nCopies(dimension, 0.0));
		this.dimension = dimension;
	}

	public Vector(double[] array, int dimension) {
		this.data = new ArrayList<>(dimension); 
		this.dimension = dimension;

		for(int i = 0; i < array.length; i++) {
			this.data.add(array[i]);
		}
	}

	public Vector(Vector newVector) {
		this.data = new ArrayList<>(newVector.getVector());
		this.dimension =  data.size();
	}

	public Vector scale(double scalar) {

		for(int i = 0; i < data.size(); i++) {
			data.set(i,data.get(i) * scalar);
		}
		return this;
	}

	public Vector add(Vector addend) {
		if(data.size() != addend.getSize()) {
			return null;
		} 


		for(int i = 0; i < data.size(); i++) {
			data.set(i,data.get(i) + addend.getDataAtIndex(i));
		}

		return this;
	}

	public static Vector Gauss_Jordan(List<Vector> vectors, int dimension, Vector constants) {
		List<Vector> transposedMatrix = Vector.transpose(vectors, dimension);

		if (transposedMatrix.size() != transposedMatrix.get(0).getSize())
			return null;
		
		Vector.rowEchelon(transposedMatrix, vectors.size(), constants);
		Vector.reducedRowEchelon(transposedMatrix, vectors.size(), constants);
		
		ArrayList<Double> zero = new ArrayList<>();
		
		for(int i=0; i<vectors.get(0).getSize(); i++)
			zero.add(0.0);
		
	
		for (int i=0; i<transposedMatrix.size(); i++) {
			if (zero.equals(transposedMatrix.get(i).getVector())) {
				return null;
			}
				
		}
			
		return constants;
	}

	public static void rowEchelon(List<Vector> vectors, int dimension, Vector constants){
		int startingRow = 0;
		for (int col = 0; col < dimension; col++) {
			int max = startingRow;
			
			while(vectors.get(max).getDataAtIndex(col) == 0 && max < vectors.size()-1)
				max++;
			
			if(vectors.get(max).getDataAtIndex(col) != 0) {
				if(max != startingRow) {
					Vector v = vectors.get(max);
					vectors.set(max, vectors.get(startingRow));
					vectors.set(startingRow, v);
					
					double temp = constants.getDataAtIndex(max);
					constants.setValue(max, constants.getDataAtIndex(startingRow));
					constants.setValue(startingRow, temp);
				}

				Vector base = new Vector(vectors.get(max));
				for(int row = startingRow+1; row < vectors.size(); row++) {
					if(vectors.get(row).getDataAtIndex(col) != 0) {
						double factor = -1 * vectors.get(row).getDataAtIndex(col);

						double first = constants.getDataAtIndex(row)*base.getDataAtIndex(col),
							   second = constants.getDataAtIndex(startingRow)*factor;

						vectors.get(row).scale(base.getDataAtIndex(col)).add(base.scale(factor));
					
						base.scale(1/factor);
						constants.setValue(row,first+second);
					}

				}

				startingRow++;
			}


		}
	}

	public static void reducedRowEchelon(List<Vector> vectors, int dimension, Vector constants){
		for(int i = vectors.size() - 1; i >= 0 ; i--) {
			int nonzeroIndex = -1;
			for(int j = 0; j < dimension; j++) {
				if(vectors.get(i).getDataAtIndex(j) != 0) {
					nonzeroIndex = j;
					break;
				}
			}

			if(nonzeroIndex != -1) {
				double factor = vectors.get(i).getDataAtIndex(nonzeroIndex);
				vectors.get(i).scale(1/factor);
				constants.setValue(i, constants.getDataAtIndex(i) / factor);

				if(i > 0) {
					for(int k = i - 1; k >= 0; k--) {
						if(vectors.get(k).getDataAtIndex(nonzeroIndex) != 0) {
							Vector base = new Vector(vectors.get(i));
							double factor2 = -1 * vectors.get(k).getDataAtIndex(nonzeroIndex);

							double first = constants.getDataAtIndex(k),
								   second = constants.getDataAtIndex(i)*factor2;

							vectors.get(k).scale(base.getDataAtIndex(nonzeroIndex)).add(base.scale(factor2));
							base.scale(1/factor2);
							constants.setValue(k,first+second);	
						}
					}
				}
			}
		}

	}

	public static int span(List<Vector> vectors, int dimension) {
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

	public static List<Vector> transpose(List<Vector> vectors, int dimension) {
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
	
	public ArrayList<Double> getVector() {
		return data;
	}

	public int getSize() {
		return data.size();
	}

	public void setValue(int index, double newValue) {
		data.set(index, newValue);
	}

	public Double getDataAtIndex(int i) {
		return data.get(i);
	}
	
	public ArrayList<Double> getData(){
		return data;
	}
	
	static int min(int a, int b) { 
		return a < b ? a : b; 

	}
}