package cn.hnu.beans;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.Vector;

import cn.hnu.actions.FileOperation;
import cn.hnu.actions.GoAnnotation;
import cn.hnu.actions.PPNHandler;



public class Network {
	
	public static int ProteinNUM ; //�����������еĵ�������Ŀ 4930
	
	public static int PPIsNUM;
	
	public static ArrayList<Node> proteins;  //�����еĵ�����
	
//	public static HashMap<Integer,Node> idx_node;  ////
	public static HashMap<String,Integer> pro_idx; //�����ʵ����id
	
//	public static HashMap<String,Node> pro_node;
	
	public static int graph[][]; //ͼ�ڽӾ���
	
	public static HashMap<String, HashSet<Integer>> goitem_pros;//go items��Ӧ�ĵ����ʼ���
	
	
	public static int GoitemsNUM ; // GO items����
	
	public static double[][] goitemssim; // GO items֮��Ĺ���������
	
	public static HashMap<String, Integer> goitems_idx; //GO itemsӳ���Ӧ��id
	
	
	public static double pro_gosim[][]; //������֮���GO������
	
}
