package cn.hnu.actions;

import java.util.List;
import java.util.TreeSet;
import java.util.Vector;


import cn.hnu.beans.Network;

public class CFOCM {
	
	
	public static void main(String[] args) throws Exception
	{
		
		///
		PPNHandler.init_ppn("ppn/dip.txt");
		/// map proteins to its corresponding GO items, and map GO items to its annotated proteins
		GoAnnotation.Pros_GoAnnotaion("go/go_note.txt");
		
		//calculate functional interdependence between each GO item pair,run only once,  
		//after that, just read matrix from file gosim.txt
		// GoAnnotation.Compute_GoItems_Similar();
		
		Network.goitemssim = FileOperation.ReadMatrixFromFile("gosim.txt");

		//generate maximum cliques
		List<Vector<TreeSet<String>>> clqs = MaxCliquesALG.GetMaxCliques(true, false);
	

		// neighbor affinity score threshold
		double threshold = 0.4;
	
		Vector<TreeSet<String>> cores = CFOCM_CORE.Detect_Complex_Cores(clqs, threshold);
			// save cores
		String prefix = "cores/un_dip_cores_";
		String suffix = ".txt";
		String corefile = prefix + String.valueOf(threshold) + suffix;
		FileOperation.SaveComplexs2File(corefile, cores);
		System.out.println("Cores saved to Path:" + corefile);

		// add attachments Append\_Core\_Attachments
		prefix = "complexes/un_dip_complexes_";
		String complexesfile =  prefix + String.valueOf(threshold) + suffix;
		CFOCM_ATTACH.Append_Core_Attachments(corefile, complexesfile);

		//evaluation
		Evaluation e = new Evaluation(complexesfile, "benchmark/mips.txt");
		e.outputFile = "result/un_dip_mips_"+ String.valueOf(threshold) + suffix;
		e.calculate();
		e = new Evaluation(complexesfile, "benchmark/sgd.txt");
		e.outputFile = "result/un_dip_sgd_"+ String.valueOf(threshold) + suffix;
		e.calculate();
		e = new Evaluation(complexesfile, "benchmark/cyc2008.txt");
		e.outputFile = "result/un_dip_408_"+ String.valueOf(threshold) + suffix;
		e.calculate();
		
	}
}
