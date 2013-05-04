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
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleInsets;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
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
	
	public static void drawHeatMap(DoubleMatrix2D matrix, String pngfile)
	{
		 double minvalue=0;
		 double maxvalue=0;
		 for (int i = 0; i < matrix.rows(); i++) {
			for (int j = 0; j < matrix.columns(); j++) {
				double temp=matrix.getQuick(i, j);
				if(!Double.isNaN(temp))
				{
					if(temp<minvalue)
						minvalue=temp;
					if(temp>maxvalue)
						maxvalue=temp;
				}

			}
		}
		 String[] tokens = pngfile.split("\\.(?=[^\\.]+$)");
		 String title=tokens[0];
		 NumberAxis numberaxis1 = new NumberAxis("Sequence");
		 numberaxis1.setRange(new Range(0, matrix.rows()));
		 NumberAxis numberaxis2 = new NumberAxis("Position");
		 numberaxis2.setRange(new Range(0, matrix.columns()));
		 DefaultXYZDataset xyzdataset = new DefaultXYZDataset();
		 xyzdataset.addSeries(title, sparseMatrix(matrix));  //here only the non-zero element will be colored
		 XYBlockRenderer xyblockrenderer = new XYBlockRenderer();
//	        LookupPaintScale lookuppaintscale = new LookupPaintScale(-1D, Double.MAX_VALUE, Color.black);
//	        lookuppaintscale.add(0D, Color.blue);
//	        lookuppaintscale.add(0.5D, Color.green);
//	        lookuppaintscale.add(1D, Color.orange);
//	        lookuppaintscale.add(2D, Color.red);
	        xyblockrenderer.setPaintScale( new LookupPaintScale(minvalue,maxvalue,Color.red));
	        

//	        xyblockrenderer.setPaintScale(getPaintScale(minvalue, maxvalue));
	       
	        
	        XYPlot xyplot = new XYPlot(xyzdataset, numberaxis2,numberaxis1, xyblockrenderer);xyplot.setBackgroundPaint(Color.lightGray);
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
