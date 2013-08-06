import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.SymbolAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleInsets;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;


public class DrawUtil {
	
	
	public static  void DrawROC( HashMap<String,XYSeries>  ROCdata, String pngfile)
	{
		  XYSeriesCollection dataset = new XYSeriesCollection();
		 for(String name: ROCdata.keySet())
		 {
			 dataset.addSeries(ROCdata.get(name));
		 }
		 JFreeChart chart = ChartFactory.createXYLineChart(
	                "ROC curve", // chart title
	                "False Positive Rate", // x axis label
	                "True Positive Rate", // y axis label
	                dataset, // data
	                PlotOrientation.VERTICAL,
	                true, // include legend
	                true, // tooltips
	                false // urls
	                );
	// NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
	        chart.setBackgroundPaint(Color.white);
	// get a reference to the plot for further customisation...
	        XYPlot plot = (XYPlot) chart.getPlot();
	        plot.setBackgroundPaint(Color.lightGray);
	        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
	        plot.setDomainGridlinePaint(Color.white);
	        plot.setRangeGridlinePaint(Color.white);
	        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
	        renderer.setShapesVisible(true);
	        renderer.setShapesFilled(true);
	// change the auto tick unit selection to integer units only...
	        NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
	        rangeAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
	        
	        ChartPanel chartPanel = new ChartPanel(chart);
	        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
	      //  setContentPane(chartPanel);
	        try {
				ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 800, 600);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

	}
	
	public static double[][] sparseMatrix(DoubleMatrix2D matrix)
	{
		IntArrayList rows = new IntArrayList();
		IntArrayList cols = new IntArrayList();
		DoubleArrayList vals = new DoubleArrayList();
		matrix.getNonZeros(rows, cols, vals);
		double[][] ret=new double[3][rows.size()];
		for (int i = 0; i < rows.size(); i++) {
			ret[0][i]=rows.get(i);
			ret[1][i]=cols.get(i);
			ret[2][i]=vals.get(i);
		}
		return ret;
	}
	
	 public static Color blend (Color color1, Color color2, double ratio)
	  {
	    float r  = (float) ratio;
	    float ir = (float) 1.0 - r;

	    float rgb1[] = new float[3];
	    float rgb2[] = new float[3];    

	    color1.getColorComponents (rgb1);
	    color2.getColorComponents (rgb2);    

	    Color color = new Color (rgb1[0] * r + rgb2[0] * ir, 
	                             rgb1[1] * r + rgb2[1] * ir, 
	                             rgb1[2] * r + rgb2[2] * ir);
	    
	    return color;
	  }
	
	 
	 	public static PaintScale getPaintScale(DoubleMatrix2D matrix)
	 	{
	 	 		
	 	 		DoubleMatrix1D vecD =matrix.viewColumn(2);
	 	 		for (int i = 3; i <matrix.columns(); i++) {
	 	 			vecD=DoubleFactory1D.dense.append(vecD, matrix.viewColumn(i));
	 			}
	 	 		DoubleMatrix1D vecSortD = vecD.viewSorted();
	 	       //... Setting PaintScale ...//
	 	 		double min=vecSortD.get(0);	
	 	 		double point1=vecSortD.get(vecSortD.size()/10);
	 	 		double point2=vecSortD.get((int) (vecSortD.size()*0.9));
	 	 		double max=vecSortD.get(vecSortD.size()-1);
	 	 	
	 	 		Color color0=Color.blue;
	 	 	    Color color1= Color.white;//blend(Color.RED,Color.white,0.01);
	 	 	    Color color2=blend(Color.RED,Color.white,0.5);
	 	 	    Color color3=Color.RED;
	 	 	    
	 	    LookupPaintScale ps = new LookupPaintScale(min, Double.MAX_VALUE, color0);
	 	    int numscale=10;
	 	   
	 	    Color purle=new Color(255, 0, 255);
	 	    ps.add(Double.NEGATIVE_INFINITY,color0);
	 	    double valPoint=min;
	 	    int num_trans=numscale;
	 	    double stepsize=(point1-min)/numscale;
	 	   for (int i = 0; i < num_trans; i++) {
	 		   ps.add(valPoint=valPoint+stepsize, blend(color1,color0,((double)i)/(num_trans)));
	 	  }
	 	   
	 	   stepsize=(point2-point1)/numscale;
	 	   if(stepsize>0)
	 	   for (int i = 0; i < num_trans; i++) {
	 		   ps.add(valPoint=valPoint+stepsize, blend(color2,color1,((double)i)/(num_trans)));
	 	  }
	 	   
	 	   stepsize=(max-point2)/numscale;
	 	   for (int i = 0; i < num_trans; i++) {
	 		   ps.add(valPoint=valPoint+stepsize, blend(color3,color2,((double)i)/(num_trans)));
	 	  }
	 	   ps.add(Double.MAX_VALUE,color3);
	 	    return ps;
	 	}
	 
	public static PaintScale getPaintScale(double min, double max)
	{
	       //... Setting PaintScale ...//
        LookupPaintScale ps = new LookupPaintScale(min, Double.MAX_VALUE, Color.gray);
        int numscale=30;
        double stepsize=(max-min)/numscale;
        Color purle=new Color(255, 0, 255);
        ps.add(Double.MIN_VALUE, Color.gray);
        double valPoint=min;
        int num_trans=numscale/3;
       for (int i = 0; i < num_trans; i++) {
    	   ps.add(valPoint=valPoint+stepsize, blend(Color.blue,Color.gray,((double)i)/(num_trans)));
	  }
       for (int i = 0; i < num_trans; i++) {
    	   ps.add(valPoint=valPoint+stepsize, blend(Color.GREEN,Color.blue,((double)i)/(num_trans)));
	  }
       for (int i = 0; i < num_trans; i++) {
    	   ps.add(valPoint=valPoint+stepsize, blend(Color.RED,Color.GREEN,((double)i)/(num_trans)));
	  }
       ps.add(Double.MAX_VALUE,Color.RED);
        return ps;
	}
	
	public static void drawHeatMap(DoubleMatrix2D matrix, String pngfile)
	{
//		 double minvalue=0;
//		 double maxvalue=0;
//		 for (int i = 0; i < matrix.rows(); i++) {
//			for (int j = 0; j < matrix.columns(); j++) {
//				double temp=matrix.getQuick(i, j);
//				if(!Double.isNaN(temp))
//				{
//					if(temp<minvalue)
//						minvalue=temp;
//					if(temp>maxvalue)
//						maxvalue=temp;
//				}
//
//			}
//		}
//		 
		 
		 String[] tokens = pngfile.split("\\.(?=[^\\.]+$)");
		 String title=tokens[0];
		 NumberAxis numberaxis1 = new NumberAxis("Sequence");
		 numberaxis1.setRange(new Range(0, matrix.rows()));
		 NumberAxis numberaxis2 = new NumberAxis("Position");
		 numberaxis2.setRange(new Range(0, matrix.columns()));
		 DefaultXYZDataset xyzdataset = new DefaultXYZDataset();
		 xyzdataset.addSeries(title, sparseMatrix(matrix));  //here only the non-zero element will be colored
		 XYBlockRenderer xyblockrenderer = new XYBlockRenderer();        

	        xyblockrenderer.setPaintScale(getPaintScale(matrix));
	       
	        
	        XYPlot xyplot = new XYPlot(xyzdataset, numberaxis1,numberaxis2, xyblockrenderer);xyplot.setBackgroundPaint(Color.white);
	        xyplot.setDomainGridlinePaint(Color.white);
	        xyplot.setRangeGridlinePaint(Color.white);
	        xyplot.setForegroundAlpha(0.66F);
//	        xyplot.setAxisOffset(new RectangleInsets(5D, 5D, 5D, 5D));
	        JFreeChart jfreechart = new JFreeChart(title, xyplot);
	        
	        try {
				ChartUtilities.saveChartAsPNG(new File(pngfile), jfreechart, 2800, 2600);
				System.out.println("Draw heatmap to file: "+pngfile);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	}
	public static void drawSignalAroundPeakCurve(DoubleMatrix2D signal_matrix,String pngfile )
	{
		
		
		//////////////prepare data points/////////////////
		XYSeriesCollection dataset = new XYSeriesCollection(); 
		XYSeries posData=new XYSeries("Positive Peak");
		XYSeries negData=new XYSeries("Negative Peak");
		double[] tempdata1=new double[signal_matrix.columns()];
		double[] tempdata2=new double[signal_matrix.columns()];
		int poscount=0;
		int negcount=0;
		for (int i = 0; i < signal_matrix.rows(); i++) {

				poscount++;
				for (int j = 0; j < signal_matrix.columns(); j++) {
					tempdata1[j]+=signal_matrix.get(i, j);
				}
			}
		
		
		for (int i = 0; i < tempdata2.length; i++) {
			posData.add(i, tempdata1[i]/poscount);
		
		}
		dataset.addSeries(posData);
	
		String[] tokens = pngfile.split("\\.(?=[^\\.]+$)");
		//////////////ploting/////////////////
		 
		 JFreeChart chart = ChartFactory.createXYLineChart(
				 tokens[0], // chart title
	                "Position Around Peak", // x axis label
	                "Signal Count", // y axis label
	                dataset, // data
	                PlotOrientation.VERTICAL,
	                true, // include legend
	                true, // tooltips
	                false // urls
	                );
	// NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
	        chart.setBackgroundPaint(Color.white);
	// get a reference to the plot for further customisation...
	        XYPlot plot = (XYPlot) chart.getPlot();
	        plot.setBackgroundPaint(Color.lightGray);
	        plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
	        plot.setDomainGridlinePaint(Color.white);
	        plot.setRangeGridlinePaint(Color.white);
	        XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
	        renderer.setShapesVisible(true);
	        renderer.setShapesFilled(true);
	// change the auto tick unit selection to integer units only...
	        NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
	        rangeAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
	        
	        ChartPanel chartPanel = new ChartPanel(chart);
	        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));
	      //  setContentPane(chartPanel);
	        try {
				ChartUtilities.saveChartAsPNG(new File(pngfile), chart, 800, 600);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	}
}
