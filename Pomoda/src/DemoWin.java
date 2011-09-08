import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.RenderingHints;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.LineBorder;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.gui.DNAStyle;
import org.biojava.bio.gui.DistributionLogo;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleInsets;




public class DemoWin {
ArrayList<PWM>  EEMPWMs=new ArrayList<PWM>(); 
ArrayList<PWM>  REMPWMs=new ArrayList<PWM>(); 
ArrayList<ArrayList<Double>> PositionDistList=new ArrayList<ArrayList<Double>>();
ArrayList<ArrayList<Double>> RankDistList=new ArrayList<ArrayList<Double>>();

public double maxdifference(ArrayList<Double> dist1,ArrayList<Double> dist2)
{
	double maxdiff=0;
	int minsize=Math.min(dist1.size(), dist2.size());
	
	for (int i = 0; i < minsize; i++) {
		double temp=Math.abs(dist1.get(i)-dist2.get(i));
		if(temp>maxdiff)
			maxdiff=temp;
	}
	return maxdiff*minsize;
}
double changeCutoff=0.05;
double changeCutoff2=0.1;

public void addPosDist(ArrayList<Double> dist)
{
	ArrayList<Double> temp=new ArrayList<Double>(dist);
	if(PositionDistList.size()==0||maxdifference(temp,PositionDistList.get(PositionDistList.size()-1))>changeCutoff2)
		PositionDistList.add(temp);
}

public void addRankDist(ArrayList<Double> dist)
{
	ArrayList<Double> temp=new ArrayList<Double>(dist);
	if(RankDistList.size()==0||maxdifference(temp,RankDistList.get(RankDistList.size()-1))>changeCutoff2)
	RankDistList.add(temp);
}
	public void addEEMmotif(PWM wm)
	{
		PWM temp=wm.Clone();
		boolean change=false;
		if(EEMPWMs.size()==0)
			change=true;
		else
		{
			if(wm.core_motiflen!=EEMPWMs.get(EEMPWMs.size()-1).core_motiflen)
				change=true;
			else if(common.PWM_Divergence(temp, EEMPWMs.get(EEMPWMs.size()-1))>changeCutoff)
				change=true;
		}
			
			
		if(change)	
		EEMPWMs.add(temp);
	}
	
	public void addREMmotif(PWM wm)
	{
		boolean change=false;
		PWM temp=wm.Clone();
		if(REMPWMs.size()==0)
			change=true;
		else
		{
			if(wm.core_motiflen!=REMPWMs.get(REMPWMs.size()-1).core_motiflen)
				change=true;
			else if(common.PWM_Divergence(temp, REMPWMs.get(REMPWMs.size()-1))>changeCutoff)
				change=true;
		}
		if(change)	
		REMPWMs.add(temp);
	}
	
	public XYPlot Array2XYPlot(ArrayList<Double> dist)
	{
		 XYSeriesCollection dataset = new XYSeriesCollection();
		  XYSeries series1 = new XYSeries("");
		//  int binsize=SearchEngine2.TotalLen/SearchEngine2.getSeqNum()/dist.size();
		  for (int i = 0; i <dist.size(); i++) {
		//	double x=binsize*i+binsize/2;//-dist.size()*this.resolution/2
		//	series1.add(x, dist.get(i));
				series1.add(i, dist.get(i));
		}
			 dataset.addSeries(series1);
			 XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
			 XYPlot plot =new XYPlot(dataset, new NumberAxis(""), null, renderer1);	 
			 XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
		        renderer.setSeriesStroke( 0, new BasicStroke(
		                4.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND
		            ));
		        renderer.setSeriesPaint(0,Color.red);
		        renderer.setShapesVisible(false);
		        renderer.setShapesFilled(true);
		        plot.getDomainAxis().setVisible(false);
		        plot.setRangeGridlinesVisible(false);
		        plot.setDomainGridlinesVisible(false);
		        plot.setOutlinePaint(Color.black);
		        plot.setOutlineVisible(true);
	        return plot;
	        
	}
	
	
	
	public JPanel PWM2Panel(PWM wm)
	{
	    int width=20;
		 int height=60;
		 JPanel wmv = new JPanel();
		 wmv.setBackground(Color.white);
		 RenderingHints hints = new  RenderingHints(RenderingHints.KEY_ANTIALIASING,  RenderingHints.VALUE_ANTIALIAS_ON);
		 wmv.setLayout(new GridLayout(1,wm.columns()));
		 for (int pos = wm.head; pos < wm.columns()-wm.tail; ++pos) {
			 DistributionLogo dl = new DistributionLogo();
			 Distribution dist = wm.getColumn(pos);
		
             dl.setRenderingHints(hints);
            // dl.setBackground(Color.white);
             
             dl.setOpaque(true);
             try {
				dl.setDistribution(dist);
			} catch (IllegalAlphabetException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
             dl.setPreferredSize(new Dimension(width, height));
             dl.setLogoPainter(new gapTextLogoPainter());
             dl.setStyle(new DNAStyle());
             dl.setBounds(pos*width, 0, width, height);
             wmv.add(dl);
		 }
		 return wmv;
	}
	public void showSEEDs(List<PWM>  seedPWMs)
	{
        JFrame frame = new JFrame("SeedFrame");
        
        frame.setLayout(new GridLayout(seedPWMs.size(),1,20,0));
        
        for (int i = 0; i < seedPWMs.size(); i++) {  
        	JPanel pan=PWM2Panel(seedPWMs.get(i));
        	pan.setBorder(new LineBorder(Color.white, 2));
        	 frame.getContentPane().add(pan);
		} 	        
        	         frame.pack();
        	         
        	         frame.setVisible(true);
	}
	
	
	public void showEEM()
	{
        JFrame frame = new JFrame("Extending EM Process");
        GridLayout layout=new GridLayout(EEMPWMs.size(),1,10,20);

        
        frame.setLayout(layout);
        
        for (int i = 0; i < EEMPWMs.size(); i++) {   	
        	JPanel pan=PWM2Panel(EEMPWMs.get(i));
        	 pan.setBorder(new LineBorder(Color.white,4));
        	 frame.getContentPane().add(pan);
		} 	        
        	         frame.pack();
        	         
        	         frame.setVisible(true);
        	         
        EEMPWMs.clear(); 	
	}
	
	
	public void showPostionDist()
	{
		 int shownum=5;
	        int remind=PositionDistList.size()/shownum;
	        if(remind<1)
	        	remind=1;
		CombinedRangeXYPlot plot = new CombinedRangeXYPlot(new NumberAxis("Probability"));
		plot.setBackgroundPaint(Color.white);
		for (int i = 0; i < PositionDistList.size(); i++) {
			if(i%remind==0)
			plot.add(Array2XYPlot(PositionDistList.get(i)), 1);
		}
		plot.setGap(50);
		NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
        rangeAxis.setLabelFont(new java.awt.Font("SansSerif", java.awt.Font.BOLD, 28));
	
		JFreeChart chart=new JFreeChart(	"Motif Position Distribution",	JFreeChart.DEFAULT_TITLE_FONT, plot, false );
		chart.setBorderPaint(Color.white);
		ChartFrame frame = new ChartFrame("Position", chart);
		frame.setPreferredSize(new Dimension(2400, 200));
		frame.pack();
		frame.setVisible(true);
		PositionDistList.clear();
	}
	
	
	public void showRankDist()
	{
		 int shownum=5;
	        int remind=RankDistList.size()/shownum;
	        if(remind<1)
	        	remind=1;
		CombinedRangeXYPlot plot = new CombinedRangeXYPlot(new NumberAxis("Probability"));
		for (int i = 0; i < 5; i++) {
			if(i%remind==0)
			plot.add(Array2XYPlot(RankDistList.get(i)), 1);
		}
		plot.setGap(50);
		plot.setBackgroundImageAlpha((float) 1.0);
		plot.setBackgroundPaint(Color.white);
		NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
	        rangeAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
	        rangeAxis.setLabelFont(new java.awt.Font("SansSerif", java.awt.Font.BOLD, 28));
	        
		JFreeChart chart=new JFreeChart(	"Motif Sequence Rank Distribution",	JFreeChart.DEFAULT_TITLE_FONT, plot, false );
		chart.setBorderPaint(Color.white);
		ChartFrame frame = new ChartFrame("Sequence Rank", chart);
		frame.setBackground(Color.white);
		frame.setPreferredSize(new Dimension(2400, 200));
		frame.pack();
		frame.setVisible(true);
		RankDistList.clear();
	}
	
	public void showREM()
	{
		  int shownum=10;
		  if(REMPWMs.size()<shownum)
			  shownum=REMPWMs.size()-1;
		  
        JFrame frame = new JFrame("Resampling EM Process");
        GridLayout layout=new GridLayout(shownum+1,1,10,20);
        
        
        frame.setLayout(layout);
       
      
        int remind=REMPWMs.size()/shownum;
        if(remind<1)
        	remind=1;
        for (int i = 0; i < REMPWMs.size(); i++) {   	
        	if(i%remind==0)
        	{
        		JPanel pan = PWM2Panel(REMPWMs.get(i));
        		  pan.setBorder(new LineBorder(Color.white,4));
        	 frame.getContentPane().add(pan);
        	}
		} 	        
        	         frame.pack();
        	         
        	         frame.setVisible(true);
        	         
       REMPWMs.clear(); 	         
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
