using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.IO;



namespace Verification
{
   public class PWM
    {
        int width;
        string name;

        public string Name
        {
            get { return name; }
            set { name = value; }
        }

        double[,] matrix;

        public int Width
        {
            get { return width; }
            set { width = value; }
        }
        public PWM(int w)
        {
            width = w;
            matrix=new double[w,4];


        }


        public static double PWM_Divergence(PWM p1, PWM p2)
        {
            double divg = 10;
            float temp;
            int overlap,overlap2;
            int aln=Alignment(p1, p2, out temp, out overlap);
            divg = temp;
            aln = Alignment(p1.RC(), p2, out temp, out overlap2);
            Console.WriteLine(aln);
            if (temp < divg)
                divg = temp;
            overlap = overlap2;
            return divg;
        }

             static int Alignment(PWM refmotif, PWM newmotif, out float alnScore, out int bestoverlap)
             {
		            float bestscore=10;
		            int bestaln=0;
                    bestoverlap = 0;
		            int size=newmotif.width;
		            int size2=refmotif.width;
                     alnScore=10;
		            int minsize=Math.Min(refmotif.width,newmotif.width);
		            for(int i=0-size+1;i<size2;i++)
		            {
			            float  curScore=0;
                        int pp = 0;
                        if (i > 0)
                            pp = i;
                        int j;
			            for(j=pp;j-pp<minsize;j++)
			            {
				            if(j>=size2||(j-i)>=size)
					            break;
				            if(j>=0)
				            {
                                double sump = 0;
                                for (int p = 0; p < 4; p++)
                                {
                                    double temp=refmotif.matrix[j,p]-newmotif.matrix[j-i,p];
                                    temp*=(temp);
                                    sump += (float)temp;
                                }
                                curScore +=(float) Math.Sqrt( sump);
				            }
			            }
                        //size2=minsize
			            int overlap=j;
                        ////i>0
                        if (i > 0)
                            overlap = j - i;
                        //if (size == minsize)
                        //{
                        //    overlap = minsize;
                        //    if (minsize + i > size2)
                        //        overlap = size2 - i;
                        //}

                        //if (i < 0)
                        //{
                        //    overlap = minsize + i;
                        //}
                        curScore /= overlap * (float)Math.Sqrt(2);
			            if(bestscore>curScore&&(overlap>=7||overlap>=minsize-2))
			            {
				            bestscore=curScore;
				            bestaln=i;
                            bestoverlap = overlap;
			            }
		            }
		            alnScore=bestscore;
		            return bestaln;
             }

        public PWM RC()
        {
            PWM rc = new PWM(width);
            rc.name = name;
            for (int i = 0; i < width; i++)
            {
                rc.matrix[width - i - 1, 0] = matrix[i, 3];
                rc.matrix[width - i - 1, 1] = matrix[i, 2];
                rc.matrix[width - i - 1, 2] = matrix[i, 1];
                rc.matrix[width - i - 1, 3] = matrix[i, 0];
            }

            return rc;
        }
       public  static List<PWM> LoadPWMFromFile(string file, int topk)
        {
            List<PWM> retlist = new List<PWM>();
            if (file.IndexOf(".pomoda") != -1)
               retlist= PomodaHandler(file);
            else if (file.IndexOf(".traw") != -1)
                retlist = TrawlerHandler(file);
            else if (file.IndexOf(".trans") != -1)
                retlist = TransfacHandler(file);
            else if (file.IndexOf(".admout") != -1)
                retlist = AmadeusHandler(file);
            else if (file.IndexOf(".wee") != -1)
                retlist = WeederHandler(file);

           if(retlist.Count>topk)
           retlist.RemoveRange(topk,retlist.Count-topk);
           return retlist;
           
        }

       static List<PWM> PomodaHandler(string file)
       {
           List<PWM> retlist = new List<PWM>();
           string line = "";
           StreamReader sr = new StreamReader(file);
           while ((line = sr.ReadLine()) != null)
           {
               if (line.IndexOf("DE\tMotif") != -1)
               {
                   string nname = line.Substring(3);
                   sr.ReadLine();
                   List<List<float>> temlist = new List<List<float>>();
                   while ((line = sr.ReadLine()) != null)
                   {
                       string[] comps=line.Split(new char[] { '\t' });
                       if (comps.Length < 4)
                           break;
                       List<float> col = new List<float>();
                       for (int i = 0; i < 4; i++)
                       {
                           col.Add(float.Parse(comps[i + 1]));
                       }
                       temlist.Add(col);
                   }

                   int w = temlist.Count;
                   if (w < 5)
                       continue;
                   PWM candidate = new PWM(w);
                   candidate.Name = nname;
                   for (int i = 0; i < w; i++)
                   {
                       for (int j = 0; j < 4; j++)
                           candidate.matrix[i, j] = temlist[i][j];
                   }
                   retlist.Add(candidate);

               }
           }
           sr.Close();
                   

           return retlist;
       }

       static List<PWM> TransfacHandler(string file)
       {
           List<PWM> retlist = new List<PWM>();
           string line = "";
           StreamReader sr = new StreamReader(file);
           string nname = "";
           while ((line = sr.ReadLine()) != null)
           {
              
               if (line.IndexOf("PO") != -1)
               {
                   
                   List<List<float>> temlist = new List<List<float>>();
                   while ((line = sr.ReadLine()) != null)
                   {
                       string[] comps = line.Split(new char[] { ',' });
                       if (comps.Length < 4)
                           break;
                       List<float> col = new List<float>();
                       float sum = 0;
                       for (int i = 0; i < 4; i++)
                       {
                           float tt = float.Parse(comps[i + 1]);
                           sum += tt;
                           col.Add(tt);
                       }
                       for (int i = 0; i < 4; i++)
                           col[i] /= sum;

                       temlist.Add(col);
                   }

                   int w = temlist.Count;
                   if (w < 5)
                       continue;
                   PWM candidate = new PWM(w);
                   candidate.Name = nname;
                   for (int i = 0; i < w; i++)
                   {
                       for (int j = 0; j < 4; j++)
                           candidate.matrix[i, j] = temlist[i][j];
                   }
                   retlist.Add(candidate);

               }
               if(line!=null)
               nname = line.Substring(2);

           }
           sr.Close();


           return retlist;
       }

       static List<PWM> TrawlerHandler(string file)
       {
           List<PWM> retlist = new List<PWM>();
           string line = "";
           StreamReader sr = new StreamReader(file);
           while ((line = sr.ReadLine()) != null)
           {
              
               while (line.IndexOf(">family") != -1)
               {
                   string nname = line.Substring(1);
                   List<List<float>> temlist = new List<List<float>>();
                   while ((line = sr.ReadLine()) != null)
                   {
                       string[] comps = line.Split(new char[] { '\t' });
                       if (comps.Length < 4)
                           break;
                       List<float> col = new List<float>();
                       float sum = 0;
                       for (int i = 0; i < 4; i++)
                       {
                           float tt = float.Parse(comps[i]);
                           sum += tt;
                           col.Add(tt);
                       }
                       for (int i = 0; i < 4; i++)
                           col[i] /= sum;

                       temlist.Add(col);
                   }
                   if (line == null)
                       break;

                   int w = temlist.Count;
                   if (w < 5)
                       continue;
                   PWM candidate = new PWM(w);
                   candidate.Name = nname;
                   for (int i = 0; i < w; i++)
                   {
                       for (int j = 0; j < 4; j++)
                           candidate.matrix[i, j] = temlist[i][j];
                   }
                   retlist.Add(candidate);
               }
           }
           sr.Close();

           return retlist;
       }

       static List<PWM> AmadeusHandler(string file)
       {
           List<PWM> retlist = new List<PWM>();
           string line = "";
           StreamReader sr = new StreamReader(file);
           while ((line = sr.ReadLine()) != null)
           {
               if (line.IndexOf("PwmMotif(") != -1)
               { 
                   int pos=line.IndexOf(':');
                   string nname = line.Substring(0, pos);
                
                 
                  line= sr.ReadLine();
                  string[] comps = line.Split(new char[] { '|' });
                  int w = comps.Length - 2;
                  PWM candidate = new PWM(w);
                  candidate.Name = nname;
                  sr.ReadLine();
                   for (int i = 0; i < 4; i++)
                   {
                       line = sr.ReadLine();
                       comps = line.Split(new char[] { '|' });
                       for (int j = 0; j < w; j++)
                       {
                           candidate.matrix[j, i] = double.Parse(comps[j + 1]);
                       }
                   }
                   retlist.Add(candidate);
               }
           }
           sr.Close();

           return retlist;
       }

       static List<PWM> WeederHandler(string file)
       {
           List<PWM> retlist = new List<PWM>();
           string line = "";
           StreamReader sr = new StreamReader(file);
           while ((line = sr.ReadLine()) != null)
           {
               string nname = "Motif_" + retlist.Count.ToString();
               if (line.IndexOf("All Occurrences") != -1)
               {
                  
                   sr.ReadLine();
                   List<List<float>> temlist = new List<List<float>>();
                   while ((line = sr.ReadLine()) != null)
                   {
                       string[] comps = line.Split(new string[] { " " },StringSplitOptions.RemoveEmptyEntries);
                       if (comps.Length < 4)
                           break;
                       List<float> col = new List<float>();
                       float sum = 0;
                       for (int i = 0; i < 4; i++)
                       {
                           float tt = float.Parse(comps[i + 2]);
                           sum += tt;
                           col.Add(tt);
                       }
                       for (int i = 0; i < 4; i++)
                           col[i] /= sum;
                       temlist.Add(col);
                   }

                   int w = temlist.Count;
                   if (w < 5)
                       continue;
                   PWM candidate = new PWM(w);
                   candidate.Name = nname;
                   for (int i = 0; i < w; i++)
                   {
                       for (int j = 0; j < 4; j++)
                           candidate.matrix[i, j] = temlist[i][j];
                   }
                   retlist.Add(candidate);

               }
           }
           sr.Close();

           return retlist;
       }



       
    }
}
