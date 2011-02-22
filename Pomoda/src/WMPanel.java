import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.RenderingHints;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.biojava.bio.*;
import org.biojava.bio.dist.*;
import org.biojava.bio.dp.*;
import org.biojava.bio.gui.*;


public class WMPanel extends JPanel {
	 private WeightMatrix wm;
     private DistributionLogo[] logos;

     public WMPanel(WeightMatrix wm) {
         super();
         this.wm = wm;
         setBackground(Color.white);

         RenderingHints hints = new 
RenderingHints(RenderingHints.KEY_ANTIALIASING, 
RenderingHints.VALUE_ANTIALIAS_ON);

         try {
             setLayout(new GridLayout(1, wm.columns()));
             logos = new DistributionLogo[wm.columns()];
             for (int pos = 0; pos < wm.columns(); ++pos) {
                 Distribution dist = wm.getColumn(pos);
                 DistributionLogo dl = new DistributionLogo();
                 dl.setRenderingHints(hints);
                 dl.setBackground(Color.white);
                 dl.setOpaque(true);
                 dl.setDistribution(dist);
                 dl.setPreferredSize(new Dimension(40, 50));
                 dl.setLogoPainter(new TextLogoPainter());
                 dl.setStyle(new DNAStyle());
                 add(dl);
                 logos[pos] = dl;
             }
         } catch (BioException ex) {
             throw new BioError(ex);
         }
     }

     public static void wmViewer(WeightMatrix wm, String message)
     {
         WMPanel wmv = new WMPanel(wm);
         JFrame frame = new JFrame("Weight matrix viewer" + ((message == 
null) ? "" : (" (" + message + ")")));
         frame.getContentPane().add(wmv);
         frame.pack();
         frame.setVisible(true);
     }
}
