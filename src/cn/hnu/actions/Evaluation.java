package cn.hnu.actions;

import java.io.*;
import java.util.*;

@SuppressWarnings("unchecked")
public class Evaluation {

	String predictedFile;
	String goldFile;

	String outputFile;

	double overlapThreshold;

	LinkedList<Vector<String>> complexList;
	LinkedList<Vector<String>> goldList;

	int n;
	int np;
	int nb;
	double Precision;
	double Recall;
	double F;
	double Accuracy;
	double Ppv;
	double Sensitivity;

	public Evaluation(String predictedFile, String goldFile) {
		this.predictedFile = predictedFile;
		this.goldFile = goldFile;
		String output = new String();
		output = predictedFile.substring(0, predictedFile.lastIndexOf("."));
		outputFile = output + "_output.txt";
		overlapThreshold = 0.2;
	}

	public void calculate() {
		// ChangeFormat(predictedFile);
		complexList = readComplexsFrom(predictedFile);
	
		System.out.println("HEEEEEEE"+complexList.size());
		goldList = readComplexsFrom(goldFile);
		
//		for(int i=0;i<goldList.size();i++)
//		{
//			Vector<String> complex = goldList.get(i);
//			for(int j=0;j<complex.size();j++)
//				System.out.print(complex.get(j)+" ");
//			System.out.println();
//		}
//		
		for(int i=0;i<complexList.size();i++)
			if(complexList.get(i).size()<3)
			{
				complexList.remove(i);
				i--;
				System.out.println("yes");
			}
		System.out.println(goldList.size()+" "+complexList.size());
		computeF(complexList, goldList);
		computeACC(complexList, goldList);
		printResult(complexList, goldList);
		print();

	}

	private void ChangeFormat(String filename) {
		BufferedReader br;
		try {
			System.out.print(filename);
			String filename2 = filename + "2";
			br = new BufferedReader(new FileReader(filename));
			String line;
			line = br.readLine();
//			String p[] = null;
			PrintWriter pw = new PrintWriter(filename2);
			while (line != null) {
				if(line.indexOf(':', 0)!=-1)
					line  =line.substring(line.indexOf(':', 0)+2);	
				pw.println(line.toUpperCase());
				line = br.readLine();
			}
			br.close();
			pw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	
	
	private LinkedList<Vector<String>> readComplexsFrom(String filename) {
		LinkedList<Vector<String>> complexs = new LinkedList<Vector<String>>();
		BufferedReader br;
		try {
			System.out.print(filename);
			br = new BufferedReader(new FileReader(filename));
			Vector<String> complex = new Vector<String>();
			String line;
			line = br.readLine();
			String p[] = null;
			while (line != null) {
//				p = split(line,'\t');
//				if(filename.contains("mips"))
					p= line.split(" |\t");
//				else
//					p= line.split(" ");
//				if(p.length>=3)
//				{
					for (int i = 0; i < p.length; i++) {
						if(p[i].trim().length()!=0)
						 complex.add(p[i]);
					}
					complexs.add(complex);
//				}
//				else
//				{
//					System.out.println("NO");
//				}
				complex = new Vector<String>();
				line = br.readLine();
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return complexs;
	}

	public String[] split(String s, char c) {
		Vector vector = new Vector(0, 1);
		String s1;
		for (s1 = s; s1.indexOf(c) != -1; s1 = s1.substring(s1.indexOf(c) + 1,
				s1.length()))
			if (s1.substring(0, s1.indexOf(c)).length() > 0)
				vector.addElement(new String(s1.substring(0, s1.indexOf(c))));

		vector.addElement(new String(s1));
		String as[] = new String[vector.capacity()];
		for (int i = 0; i < vector.capacity(); i++)
			as[i] = (String) vector.elementAt(i);

		return as;
	}


	public void computeF(LinkedList<Vector<String>> complex,
			LinkedList<Vector<String>> gold) {
		double score = 0.0D;
		for (int i = 0; i < complex.size(); i++) {
			for (int j = 0; j < gold.size(); j++) {
				score = getSimilarity(complex.get(i), gold.get(j));
				if (score >= overlapThreshold) {
					np++;
					break;
				}
			}
		}

		for (int i = 0; i < gold.size(); i++) {
			for (int j = 0; j < complex.size(); j++) {
				score = getSimilarity(complex.get(j), gold.get(i));

				if (score >= overlapThreshold) {
					nb++;
					break;
				}
			}
		}
		Precision = Double.valueOf(np) / complex.size();
		Recall = Double.valueOf(nb) / gold.size();
		F = 2 * Precision * Recall / (Precision + Recall);
	}

	public void computeACC(LinkedList<Vector<String>> complex,
			LinkedList<Vector<String>> gold) {
		int maxnum = 0;
		int sum = 0;
		int goldprotein = 0;
		int proteins = 0;
		int num = 0;
		for (int j = 0; j < gold.size(); j++) {
			goldprotein += gold.get(j).size();
		}
		for (int i = 0; i < complex.size(); i++) {
			for (int j = 0; j < gold.size(); j++) {
				proteins += getIntersectionSize(complex.get(i), gold.get(j));
			}
		}

		int i = 0;
		while (i < gold.size()) {
			maxnum = 0;
			for (int j = 0; j < complex.size(); j++) {
				num = getIntersectionSize(complex.get(j), gold.get(i));
				if (num > maxnum)
					maxnum = num;
			}
			sum += maxnum;
			i++;
		}
		Sensitivity = Double.valueOf(sum) / goldprotein;

		i = 0;
		sum = 0;
		while (i < complex.size()) {
			maxnum = 0;
			for (int j = 0; j < gold.size(); j++) {
				num = getIntersectionSize(complex.get(i), gold.get(j));
				if (num > maxnum)
					maxnum = num;
			}
			sum += maxnum;
			i++;
		}
		Ppv = Double.valueOf(sum) / proteins;
		Accuracy = Math.sqrt(Ppv * Sensitivity);
	}

	public void printResult(LinkedList<Vector<String>> complex,
			LinkedList<Vector<String>> gold) {
		try {
			FileOutputStream fileoutputstream = new FileOutputStream(outputFile);
			DataOutputStream dataoutputstream = new DataOutputStream(
					fileoutputstream);
			
			String ss = String.valueOf(complex.size())+" "+ String.valueOf(np)+" "+
					String.valueOf(nb)+" "+String.format("%.4f",Precision)+" "+
					String.format("%.4f",Recall)+" "+ String.format("%.4f",F)+" "+
					String.format("%.4f",Sensitivity)+" "+String.format("%.4f",Ppv)+" "+ String.format("%.4f",Accuracy)+"\n";
			dataoutputstream.writeBytes(ss);
			dataoutputstream.close();
		} catch (Exception exception) {
			System.out.println("ÊÇÂð+++++++++" + exception);
		}
	}

	public double getSimilarity(Vector set1, Vector set2) {
		double num = getIntersectionSize(set1, set2);
		if (num == 0)
			return 0;
		return num * num * 1.0 / set1.size() / set2.size();
	}

	public int getIntersectionSize(Vector set1, Vector set2) {
		int num = 0;
		Iterator iter = null;
		Vector hSet = null;
		if (set1.size() > set2.size()) {
			iter = set2.iterator();
			hSet = set1;
		} else {
			iter = set1.iterator();
			hSet = set2;
		}
		for (; iter.hasNext();) {
			if (hSet.contains(iter.next())) {
				num++;
			}
		}
		return num;
	}

	public double qSimilary(Hashtable ht1, Hashtable ht2) {
		int n = 0;
		double t = 0.0d;
		for (Enumeration enumeration1 = ht1.keys(); enumeration1
		.hasMoreElements();) {
			Integer integer1 = (Integer) enumeration1.nextElement();
			if (ht2.containsKey(integer1)) {
				n++;
			}
		}
		t = (n + 1) / Math.sqrt(ht1.size() * ht2.size());
		return t;
	}

	public void print() {
		System.out.println("complex:" + complexList.size());
		System.out.println("gold: " + goldList.size());
		System.out.println("np:" + np);
		System.out.println("nb:" + nb);
		System.out.println("Precision: " + Precision);
		System.out.println("Recall: " + Recall);
		System.out.println("F: " + F);
	    System.out.println("Sensitivity: " + Sensitivity);
		System.out.println("Ppv: " + Ppv);
		System.out.println("Accuracy: " + Accuracy);
	}

	public static void main(String args[]) {

	}

}
