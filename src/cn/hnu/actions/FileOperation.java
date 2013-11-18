package cn.hnu.actions;



import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.Vector;

public class FileOperation {

	public static Vector<TreeSet<String>> ReadComplexsFromFile(String filename) {
		Vector<TreeSet<String>> clqs = new Vector<TreeSet<String>>();
		try {
			FileReader f = new FileReader(filename);
			BufferedReader in = new BufferedReader(f);
			String s = in.readLine();
			while (s != null) {
				String[] temp = s.split(" |\t");
				TreeSet<String> clq = new TreeSet<String>();
				for (int i = 0; i < temp.length; i++) {
					clq.add(temp[i].toUpperCase());
				}
				clqs.add(clq);
				s = in.readLine();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(clqs.size());
		HashSet<TreeSet<String>> set = new HashSet<TreeSet<String>>();
		set.addAll(clqs);
		//System.out.println(set.size());
		clqs = new Vector<TreeSet<String>>();
		clqs.addAll(set);
//		System.out.println(clqs.size());
		return clqs;
		
	}
	
	public static HashSet<TreeSet<String>> ReadComplexsFromFile2(String filename) {
		HashSet<TreeSet<String>> clqs = new HashSet<TreeSet<String>>();
		try {
			FileReader f = new FileReader(filename);
			BufferedReader in = new BufferedReader(f);
			String s = in.readLine();
			while (s != null) {
				String[] temp = s.split(" |\t");
				TreeSet<String> clq = new TreeSet<String>();
				for (int i = 0; i < temp.length; i++) {
					clq.add(temp[i].toUpperCase());
				}
				clqs.add(clq);
				s = in.readLine();
			}
		} catch (Exception e) {
		}
		return clqs;
	}
	
	public static void SaveComplexs2File(String filename,
			Vector<TreeSet<String>> clqs) {
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(filename);
			for (int i = 0; i < clqs.size(); i++) {
				String[] ps =  new String[clqs.elementAt(i).size()];
				clqs.elementAt(i).toArray(ps);
				int j;
				for(j=0;j<ps.length-1;j++)
				{
					pw.print(ps[j] + " ");
				}
				pw.println(ps[ps.length-1]);
			}
			pw.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static double[][] ReadMatrixFromFile(String filename) {
		double[][] A = null;
		try {
			String line;
			boolean first = true;
			int len = 0;
			int j=0;
			BufferedReader br = new BufferedReader(new FileReader(filename));
			while((line=br.readLine())!=null)
			{
				String[] temp = line.split(" ");
				if(first)
				{
					len = temp.length;
					A = new double [len][len];
					first =false;
				}
				for(int i=0;i<temp.length;i++)
				{
					A[j][i]= Double.parseDouble(temp[i]);
				}
				j++;
			}
			br.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return A;
	}
	
	
	public static void SaveMatrix2File(double A[][],String filename) {
		
		PrintWriter pw = null;
		try {
			pw = new PrintWriter(filename);
//			System.out.println(A.length + " " + A[0].length);
			for (int i = 0; i < A.length; i++) {

				for (int j = 0; j < A[0].length-1; j++) {
					pw.print(A[i][j]+" ");
					// System.out.print(A[i][j]+" ");
				}
				pw.println(A[i][A[0].length-1]);
			}
			pw.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
}
