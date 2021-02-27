package ClassicalStat;

import ClassicalStat.*;
import java.util.*;

/**  Class for combinatory calculations 
*/
public class Combinatory
{

	private RandomGen rand;

	//region Constructor

	/**  Constructor 
	*/
	public Combinatory()
	{
		rand = new RandomGen();
	}

	//endregion

	//region Basic Combinatory

	/**  Combinations of n from g 
	 @param n total number of elements 
	 @param g group of elements 
	 @return  number of combinations 
	*/
	public final double Combinations(int n, int g)
	{
		return Factorial(n) / (Factorial(g) * Factorial(n - g));
	}

	/**  Permutations of n  
	 @param n number of elements to permutate 
	 @return  number of permutations 
	*/
	public final double Permutations(int n)
	{
		return Factorial(n);
	}

	/**  Factorial of n 
	 @param n number to calculate factorial 
	 @return  factorial value 
	*/
	public final double Factorial(int n)
	{
		if (n < 0)
		{
			throw new RuntimeException("Strings.n_must_be_greater_than_zero");
		}
		if (n == 0)
		{
			return 1;
		}
		if (n == 1)
		{
			return 1;
		}
		else
		{
			return n * Factorial(n - 1);
		}
	}

	//endregion

	//region Sums Combinatory

	//region Sums Combinations

	/**  Calculate combinations of n consecutive values from the same list with replacement that sum a particular total 
	 @param values possible ordered consecutive values of summands 
	 @param n number of summands 
	 @param k sum 
	 @param order if list of values should be previously ordered 
	 @return  list of m lists of n summands that sum df
	*/
	public final ArrayList<int[]> SumsCombinations(ArrayList<Integer> values, int n, int k, boolean order)
	{
		if (order)
		{
			Collections.sort(values);
		}

		ArrayList<int[]> sumsCombinations = new ArrayList<int[]>();
		int[] comb = new int[n];
		SumsCombinations(sumsCombinations, values, n, comb, 0, 0, k);
		return sumsCombinations;
	}

	private void SumsCombinations(ArrayList<int[]> sumsCombs, ArrayList<Integer> values, int n, int[] comb, int top, int actSum, int sum)
	{

		if (sum < values.get(0) * comb.length || sum > values.get(values.size() - 1) * comb.length)
		{
			return;
		}

		if (top == n - 1)
		{
			if (actSum + values.get(0) <= sum && actSum + values.get(values.size() - 1) >= sum)
			{
				comb[top] = sum - actSum;
				int[] clone = new int[n];
				System.arraycopy(comb, 0, clone, 0, comb.length);
				sumsCombs.add(clone);
			}
		}
		else
		{
			for (int i = 0;i < values.size();i++)
			{
				comb[top] = values.get(i);
				if ((n - top - 1) * values.get(0) <= sum - actSum - values.get(i) && (n - top - 1) * values.get(values.size() - 1) >= sum - actSum - values.get(i))
				{
					SumsCombinations(sumsCombs, values, n, comb, top + 1, actSum + values.get(i), sum);
				}
			}
		}
	}

	//endregion

	//region Sums Variations

	/**  Calculate variations of n consecutive values from the same list with replacement that sum a particular total quantity 
	 @param values possible ordered consecutive values of summands 
	 @param n number of summands 
	 @param k sum 
	 @param order if list of values should be previously ordered 
	 @return  list of m lists of n summands that sum df
	*/
	public final ArrayList<int[]> SumsVariations(ArrayList<Integer> values, int n, int k, boolean order)
	{
		if (order)
		{
			Collections.sort(values);
		}

		ArrayList<int[]> sumsVariations = new ArrayList<int[]>();
		int[] var = new int[n];
		SumsVariations(sumsVariations, values, n, var, 0, 0, k, 0);
		return sumsVariations;
	}

	private void SumsVariations(ArrayList<int[]> sumsVars, ArrayList<Integer> values, int n, int[] var, int top, int actSum, int sum, int firstIndex)
	{

		if (sum < values.get(0) * var.length || sum > values.get(values.size() - 1) * var.length || firstIndex > values.size() - 1)
		{
			return;
		}

		if (top == n - 1)
		{
			if (actSum + values.get(firstIndex) <= sum && actSum + values.get(values.size() - 1) >= sum)
			{
				int diff = sum - actSum;
				for (int i = firstIndex;i < values.size();i++)
				{
					if (values.get(i).equals(diff))
					{
						var[top] = diff;
						int[] clone = new int[n];
						System.arraycopy(var, 0, clone, 0, var.length);
						sumsVars.add(clone);
					}
				}
			}
		}
		else
		{
			for (int i = firstIndex;i < values.size();i++)
			{
				var[top] = values.get(i);
				if ((n - top - 1) * values.get(0) <= sum - actSum - values.get(i) && (n - top - 1) * values.get(values.size() - 1) >= sum - actSum - values.get(i))
				{
					SumsVariations(sumsVars, values, n, var, top + 1, actSum + values.get(i), sum, i);
				}
			}
		}
	}

	/**  Calculate variations of n consecutive values from the same list with replacement that sum less or equal a particular total quantity 
	 @param values possible ordered consecutive values of summands 
	 @param n number of summands 
	 @param k sum 
	 @param order if list of values should be previously ordered 
	 @return  list of m lists of n summands that sum df
	*/
	public final HashMap<Integer, ArrayList<int[]>> SumsVariationsUntil(ArrayList<Integer> values, int n, int k, boolean order)
	{
		if (order)
		{
			Collections.sort(values);
		}

		HashMap<Integer, ArrayList<int[]>> sumsVarsDict = new HashMap<Integer, ArrayList<int[]>>();
		int[] var = new int[n];
		SumsVariationsUntil(sumsVarsDict, values, n, var, 0, 0, k, 0);
		return sumsVarsDict;
	}

	private void SumsVariationsUntil(HashMap<Integer,ArrayList<int[]>> sumsVarsDict, ArrayList<Integer> values, int n, int[] var, int top, int actSum, int sum, int firstIndex)
	{

		if (sum < values.get(0) * var.length || sum > values.get(values.size() - 1) * var.length || firstIndex > values.size() - 1)
		{
			return;
		}

		if (top == n - 1)
		{
			if (actSum + values.get(firstIndex) <= sum)
			{
				int val = 0;
				double max = sum - actSum;
				int i = firstIndex;
				int[] clone;
				int finalSum;
				while (val < max && i < values.size())
				{
					val = values.get(i);
					finalSum = actSum + val;
					var[top] = values.get(i);
					clone = new int[n];
					System.arraycopy(var, 0, clone, 0, var.length);
					if (val < max)
					{
						if (!sumsVarsDict.containsKey(finalSum))
						{
							sumsVarsDict.put(finalSum, new ArrayList<int[]>());
						}
						sumsVarsDict.get(finalSum).add(clone);
					}
					i++;
				}
			}
		}
		else
		{
			for (int i = firstIndex;i < values.size();i++)
			{
				var[top] = values.get(i);
				if ((n - top - 1) * values.get(0) <= sum - actSum - values.get(i))
				{
					SumsVariationsUntil(sumsVarsDict, values, n, var, top + 1, actSum + values.get(i), sum, i);
				}
			}
		}
	}

	//endregion

	//endregion

	//region Permutations

	/**  obtain a permutation of indexes 
	 @param indexes list of indexes 
	 @return  permutated list 
	*/
	public final ArrayList<Integer> Permutation(List<Integer> indexes)
	{
		ArrayList<Integer> perm = new ArrayList<Integer>();
		while (indexes.size() > 0)
		{
			int i = rand.NextInt(0, indexes.size() - 1);
			perm.add(indexes.get(i));
			indexes.remove(i);
		}
		return perm;
	}

	/**  obtain a permutation of indexes of an interval 
	 @param indexes list of indexes 
	 @param min min index 
	 @param max max index 
	 @return  permutated list 
	*/
	public final ArrayList<Integer> Permutation(int min, int max)
	{
		ArrayList<Integer> indexes = new ArrayList<Integer>();
		for (int i = min; i <= max; i++)
		{
			indexes.add(i);
		}
		return Permutation(indexes);
	}

	public final ArrayList<Integer> PermutationNoRepeat(List<Integer> indexes)
	{
		ArrayList<Integer> perm = new ArrayList<Integer>();
		int selIndex;
		while (indexes.size() > 0)
		{
			selIndex = rand.NextInt(0, indexes.size() - 1);
			perm.add(indexes.get(selIndex));
			indexes.remove(selIndex);
		}
		return perm;
	}

	public final ArrayList<Integer> PermutationNoRepeat(int size)
	{
		ArrayList<Integer> indexes = new ArrayList<Integer>();
		for (int i = 0; i < size; i++)
		{
			indexes.add(i);
		}
		return PermutationNoRepeat(indexes);
	}

	//endregion

}