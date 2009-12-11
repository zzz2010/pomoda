using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Windows.Forms;
using System.IO;

namespace Verification
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            //getResultFromFile("");
        }
        Dictionary<string, int> SeqNameNoMap = new Dictionary<string, int>();
        private void button1_Click(object sender, EventArgs e)
        {
            listBox1.Items.Clear();
            dataGridView1.Rows.Clear();
            compareResult(textBox_Ans.Text, textBox1.Text);
        }

        void compareResult(string ansDIR, string resDIR)
        {
            string[] ansList = Directory.GetFiles(ansDIR);
            for (int i = 0; i < ansList.Length; i++)
            {
                compareResultOne(ansList[i], resDIR);
            }
        }
        private delegate void StringToListBox(String item);

        private void AddStringToListBox(String item)
        {
            //listBox1.Items.Add(item);
            Regex rg = new Regex(@"select_MAKEW(?<ml>\d+)S(?<cns>\d)A(?<adn>\d+)V(?<cv>\d+)L(?<sl>\d+)N(?<sn>\d+)\S+\s+(?<spec>[\d,\.]+)\s+(?<rcall>\d+)");
            Match mm = rg.Match(item);
            ResultRow rr = new ResultRow();
            rr.MotifLength =( mm.Groups["ml"].Value);
            rr.Conservation = (mm.Groups["cns"].Value);
            rr.Abundance =int.Parse(  mm.Groups["adn"].Value);
            if (rr.Abundance == 2)
                rr.Abundance = 0.1;
            else
                rr.Abundance = 0.5;
            rr.CenterVariance=int.Parse( mm.Groups["cv"].Value);
            rr.SeqNum = int.Parse(mm.Groups["sn"].Value);
            rr.Specificity = double.Parse(mm.Groups["spec"].Value);
            rr.Sensitivity = double.Parse(mm.Groups["rcall"].Value) / (rr.SeqNum*rr.Abundance);
            //dataGridView1.Rows.Add();
            bindingSource1.Add(rr);
            dataGridView1.Refresh();

        }
        void compareResultOne(string ansfile, string resDIR)
        {
            List<int> ansDD = new List<int>();
            StreamReader sr = new StreamReader(ansfile);
            string line = "";
            while ((line = sr.ReadLine()) != null)
            {
                if (line != "")
                {
                    string[] comps=line.Split(new char[] { ',' });
                    int seq = int.Parse(comps[0]);
                    int pos = int.Parse(comps[1]);
                    ansDD.Add((seq + 1) * 10000 + pos);//
                }
            }
            sr.Close();
            int prefid=ansfile.IndexOf(".fa");
            int prefid2 = ansfile.LastIndexOf('\\');
            string prefix = ansfile.Substring(prefid2 + 1, prefid -prefid2);
            string[] predList = Directory.GetFiles(resDIR);
            for (int i = 0; i < predList.Length; i++)
            {
                if (predList[i].IndexOf(prefix) != -1)
                {
                    List<List<int>> predSet = getResultFromFile(predList[i]);
                    double mp;
                    int ppos = predList[i].LastIndexOf('\\');
                    string resname = predList[i].Substring(ppos + 1, predList[i].Length - ppos - 1);
                    int maxcnt;
                    //if (resname.IndexOf("W8S1A1V50.fasta.wee") != -1)
                    //    resname = resname;
                      
                    cout2ListMaxOverlap(ansDD, predSet,out mp, out maxcnt);

                    resname += "\t"+mp.ToString()+"\t"+maxcnt.ToString();

                    //if (resname.IndexOf(".bp") != -1)
                    //listBox1.BeginInvoke(
                    //new StringToListBox(AddStringToListBox),
                    // new object[] { resname });
                    dataGridView1.BeginInvoke(new StringToListBox(AddStringToListBox),
                     new object[] { resname });
                }
            }
            

        }

        void cout2ListMaxOverlap(List<int> ans, List<List<int>> pred,out double maxprecision,out int maxcnt)
        {
            maxprecision = 0;
            maxcnt = 0;
            int cnt=0;
            ans.Sort();

            for (int kk = 0; kk < pred.Count; kk++)
            {
                List<int> plist2 = pred[kk];
                plist2.Sort();
                List<int> plist1 = ans;
                int SEQLEN = 10000;
                int i, j;
                i = j = 0;
                double countIntersect = 0;
                double countUnion = 0;
                while (i < plist1.Count && j < plist2.Count)
                {

                    int seqnum1 = plist1[i] / SEQLEN;
                    int pos1 = plist1[i] % SEQLEN;
                    int seqnum2 = plist2[j] / SEQLEN;
                    int pos2 = plist2[j] % SEQLEN;
                    if (seqnum1 == seqnum2 && Math.Abs((int)(pos2 - pos1)) < 10 / 2)
                    {
                        countIntersect += 1;
                        countUnion -= 1;
                    }
                    if (plist1[i] < plist2[j])
                        i++;
                    else
                        j++;
                }

                double tempPrecision = countIntersect / plist2.Count;
                if (tempPrecision > maxprecision)
                {
                    maxprecision = tempPrecision;

                    maxcnt = (int)countIntersect;
                }
            }

            return;
        }

        List<List<int>> getResultFromFile(string filename)
        {
            Dictionary<string, Regex> formatTable = new Dictionary<string, Regex>();
            formatTable["out"] = new Regex(@">SEQ_(?<SEQ>\d+)\s(?<RC>[b,f])(?<POS>\d+)\s(?<SITE>[A,C,G,T]+)");
            formatTable["mds_out"] = new Regex(@">(?<SEQ>\S+)\s\w+\s[\w+\s]*(?<RC>[b,f])(?<POS>\d+)\s(?<SITE>[A,C,G,T]+)");
            formatTable["bpout"] = new Regex(@">SEQ_(?<SEQ>\d+)\s\w+\s[\w+\s]*(?<RC>[r,f])(?<POS>\d+)\s\S+\s(?<SITE>[A,C,G,T]+)");
            formatTable["annspec_out"] = new Regex(@"\w+\s+\S+\s+\S+\s+(?<POS>\d+)\s+(?<SITE>[A,C,G,T]+)\s+\S+\s+>SEQ_(?<SEQ>\d+)");
            formatTable["wee"] = new Regex(@"(?<SEQ>\d+)\s+(?<RC>[+,-])\s+(?<SITE>[\[,\],A,C,G,T]+)\s+(?<POS>\d+)\s+\(");
            formatTable["txt"] = new Regex(@"(?<SEQ>\S+)\s+(?<POS>\d+)\s+\d+\.\d+e");
            formatTable["meme_out"] = new Regex(@"(?<SEQ>\S+)\s+(?<RC>[+,-])\s+(?<POS>\d+)\s+\d+\.\d+e");
            //formatTable["wee"] = new Regex(@"(?<SEQ>\d+)\s+(?<RC>[+,-])\s+(?<SITE>[A,C,G,T]+)\s+(?<POS>\d+)\s+\(");

            Regex t = formatTable["txt"];
            Match mm = t.Match("MACS_LNCAPARDHT_001776      179  5.25e-06 TCAACTTATT CTCCTGCC TCAGTTTTCC");
            string ss = mm.Groups["SEQ"].Value;
            string ps = mm.Groups["POS"].Value;
            //string bs = mm.Groups["SITE"].Value;
            //string rc = mm.Groups["RC"].Value;
            List<List<int>> ret = new List<List<int>>();
            StreamReader sr = new StreamReader(filename);

            string line = "";
            int pos = filename.LastIndexOf('.');
            string sufix=filename.Substring(pos + 1, filename.Length - pos-1);
            //if (filename.IndexOf("W8S1A1V50.fasta.wee") != -1)
            //    sufix = sufix;
            Regex matcher = formatTable[sufix];
            bool newList = true;
            int listCnt = -1;
            while ((line = sr.ReadLine()) != null)
            {
                if (matcher.IsMatch(line))
                {
                    if (newList)
                    {
                        listCnt++;
                        ret.Add(new List<int>());
                        newList = false;
                    }
                    if(!newList)
                    {
                        //parse the line
                        //string[] comps=line.Split(new char[] {' ','\t',',' });
                        //int pos= line.IndexOf("SEQ_") + 4;
                        Match mat = matcher.Match(line);
                        Position pp = new Position();
                       
                        pp.POS = int.Parse(mat.Groups["POS"].Value);
                        int seqqq;
                        string seqstr = mat.Groups["SEQ"].Value.Replace("SEQ_", "");
                        int cut = seqstr.Length;
                        if (cut > 25)
                            cut = 25;
                        if (int.TryParse(seqstr, out seqqq))
                            pp.SEQ = seqqq;
                        else
                        pp.SEQ = SeqNameNoMap[mat.Groups["SEQ"].Value.Substring(0,cut)];
                        if (mat.Groups["RC"].Value == "r" || mat.Groups["RC"].Value == "b" || mat.Groups["RC"].Value == "-")
                            pp.POS = (int)numericUpDown1.Value - pp.POS;
                        //if(Math.Abs(pp.POS-numericUpDown2.Value/2)<numericUpDown1.Value/2)
                        ret[listCnt].Add(pp.SEQ * 10000 + pp.POS);

                    }
                }
                else
                {
                    newList = true;
                }
            }

            sr.Close();
            return ret;
        }

        void setUpSeqMap(string filename)
        {
            SeqNameNoMap.Clear();
            StreamReader sr = new StreamReader(filename);
            string line = "";
            int count = 1;
            while ((line = sr.ReadLine()) != null)
            {
                if (line != "" && line[0] == '>')
                {
                    int cut = line.Length - 1;
                    if (cut > 24)
                        cut = 24;
                    SeqNameNoMap[line.Substring(1, cut)] = count;
                    count++;
                }
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            string[] resfiles=Directory.GetFiles(textBox1.Text);
            for (int i = 0; i < resfiles.Length; i++)
            {
                int prefid = resfiles[i].LastIndexOf(".");
                int prefid2 = resfiles[i].LastIndexOf('\\');
                string prefix = resfiles[i].Substring(prefid2 + 1, prefid - prefid2 - 1);
                //if (prefix.IndexOf("-") != -1)
                //{
                //    int predid3 = prefix.IndexOf('-');
                //    prefix = prefix.Substring(0, predid3);
                //}
                //if (prefix.IndexOf(".fasta") == -1)
                //    prefix += ".fasta";
                setUpSeqMap(textBox2.Text + "\\" + prefix);
                Directory.CreateDirectory("./"+prefix + "Pos");
                 List<List<int>> posSet=getResultFromFile(resfiles[i]);
                 string path = "./" + prefix + "Pos/";
                 for (int j = 0; j < posSet.Count; j++)
                 {
                     StreamWriter sw = new StreamWriter(path +"Motif" + j.ToString() + ".pos");
                     List<int> posList = posSet[j];
                     posList.Sort();
                     for (int k = 0; k < posList.Count; k++)
                     {
                         int pos = posList[k] % 10000;
                         int seq = posList[k] / 10000 -1;
                         pos = pos - (int)numericUpDown2.Value/2;
                         sw.WriteLine(seq.ToString() + "," + pos.ToString());
                     }
                     sw.Close();
                 }
            }
        }

        private void Form1_Load(object sender, EventArgs e)
        {

          
        }

        private void button3_Click(object sender, EventArgs e)
        {
            List<PWM> ansList = PWM.LoadPWMFromFile("anslist.trans", 100000);
            string[] predList = Directory.GetFiles(textBox1.Text);
            bindingSource1.Clear();
            for (int i = 0; i < predList.Length; i++)
            {
                int pos1 = predList[i].LastIndexOf('\\');
                string dataname = predList[i].Substring(pos1+1);

                
                List<PWM> predPWMs = PWM.LoadPWMFromFile(predList[i], (int)numericUpDown3.Value);
                 for (int j = 0; j < predPWMs.Count; j++)
                 {
                     double bestscore = 10;
                     int ansindex = 0;
                     for (int k = 0; k < ansList.Count; k++)
                     {

                         int pos2 = ansList[k].Name.IndexOf('_');
                         if (pos2 == -1)
                             pos2 = ansList[k].Name.Length;
                         if (dataname.IndexOf(ansList[k].Name.Substring(0, pos2)) == -1&&checkBox1.Checked)
                             continue;
                         //if (j == 7 && ansList[k].Name.IndexOf("SP1") != -1)
                         //    j=j;
                         double divg = PWM.PWM_Divergence(predPWMs[j], ansList[k]);

                         if (bestscore > divg)
                         {
                             bestscore = divg;
                             ansindex = k;
                         }
                          if (divg < 0.24)
                         {
                             ResultRow rr = new ResultRow();
                             rr.MotifLength = dataname + predPWMs[j].Name;
                             rr.Conservation = ansList[k].Name;
                             if (divg < 0.12)
                                 rr.Abundance = divg;
                             else if (divg < 0.18)
                                 rr.CenterVariance = divg;
                             else
                                 rr.SeqNum = divg;
                             bindingSource1.Add(rr);
                             //dataGridView1.Refresh();
                         }
                     }
                     if (bestscore < 0.24)
                     {
                         ResultRow rr = new ResultRow();
                         rr.MotifLength = dataname + predPWMs[j].Name;
                         rr.Conservation = ansList[ansindex].Name;
                         if (bestscore < 0.12)
                             rr.Abundance = bestscore;
                         else if (bestscore < 0.18)
                             rr.CenterVariance = bestscore;
                         else
                             rr.SeqNum = bestscore;
                         bindingSource1.Add(rr);
                         
                     }
                     
                 }
                 dataGridView1.Refresh();
                 
            }

        }

        private void button4_Click(object sender, EventArgs e)
        {
            bindingSource1.Clear();
            StreamReader sr = new StreamReader("log.txt");
            StreamReader sr2 = new StreamReader("AnsconserReg.txt");
            List<string> knowSeeds=new List<string>();
            Dictionary<string, bool> knowmotifs = new Dictionary<string,bool>();
            string line = "";
            while ((line = sr2.ReadLine()) != null)
            {
                knowSeeds.Add(line);
                int pp = line.IndexOf('\t');
                knowmotifs[line.Substring(pp + 1, line.Length - pp - 1)]=false;

            }

            int redundancy = 0;
            int matchingSeed = 0;
            bool maflag = false;
            int scnt = 0;
            while ((line = sr.ReadLine()) != null)
            {
                if (line.IndexOf("top:") != -1)
                {
                    scnt++;
                    if (scnt > numericUpDown3.Value)
                        break;
                    maflag = false;
                    string temp = line.Substring(line.IndexOf('\t') + 1,16);
                    temp = temp.Substring(0, temp.IndexOf(':'));
                    for (int i = 0; i < knowSeeds.Count; i++)
                    {
                        if (knowSeeds[i].IndexOf(temp) != -1)
                        {
                            ResultRow rr = new ResultRow();
                            rr.MotifLength = temp;
                            int pp=knowSeeds[i].IndexOf('\t');
                            rr.Conservation = knowSeeds[i].Substring(pp+1,knowSeeds[i].Length-pp-1);
                            if (knowmotifs.ContainsKey(rr.Conservation))
                            {
                                if (!maflag)
                                {
                                    matchingSeed++;
                                    maflag = true;
                                }
                                redundancy++;
                                knowmotifs[rr.Conservation] = true;
                            }
                            bindingSource1.Add(rr);
                            //dataGridView1.Refresh();
                        }
                    }
                }
                else
                {
                    if (line.IndexOf("windowsize") != -1)
                        break;
                }
            }

            Dictionary<string, bool>.Enumerator iter = knowmotifs.GetEnumerator();
            int cnt = 0;
            int uniquemotif=0;
            HashSet<string> family = new HashSet<string>();
            string lastfamily = "-----";
            while (iter.MoveNext())
            {
                if (iter.Current.Value)
                {
                    ResultRow rr = new ResultRow();
                    rr.MotifLength = iter.Current.Key;
                    cnt++;
                    bindingSource1.Add(rr);
                    if (!family.Contains(rr.MotifLength.Substring(0, 3)))
                    {
                        uniquemotif++;
                        family.Add(rr.MotifLength.Substring(0, 3));
                        
                    }
                    dataGridView1.Refresh();
                }
            
            }
             ResultRow r = new ResultRow();
             r.MotifLength = ((double)((double)cnt / knowmotifs.Count)).ToString();
             r.Sensitivity = (double)cnt / uniquemotif;
             r.Abundance = (double)redundancy / cnt;
             r.CenterVariance = (double)redundancy / matchingSeed;
             bindingSource1.Add(r);
             dataGridView1.Refresh();
            

            sr.Close();
            sr2.Close();

        }

        private void button5_Click(object sender, EventArgs e)
        {
            string[] comps = textBox3.Text.Split(new char[] { ' '});
            StreamReader sr = new StreamReader("anslistAll.trans");
            StreamWriter sw = new StreamWriter("anslist.trans");
             string line = "";
             while ((line = sr.ReadLine()) != null)
             { 
                
             
             }
        }
        
    }


    class Position
    {
        int sEQ;

        public int SEQ
        {
            get { return sEQ; }
            set { sEQ = value; }
        }
        int pOS;

        public int POS
        {
            get { return pOS; }
            set { pOS = value; }
        }
    }

    class ResultRow
    {
        string motifLength;

        public string MotifLength
        {
            get { return motifLength; }
            set { motifLength = value; }
        }

       string conservation;

        public string Conservation
        {
            get { return conservation; }
            set { conservation = value; }
        }

        double abundance;

        public double Abundance
        {
            get { return abundance; }
            set { abundance = value; }
        }
        double centerVariance;

        public double CenterVariance
        {
            get { return centerVariance; }
            set { centerVariance = value; }
        }
        double seqNum;

        public double SeqNum
        {
            get { return seqNum; }
            set { seqNum = value; }
        }

        
        double specificity;

        public double Specificity
        {
            get { return specificity; }
            set { specificity = value; }
        }
        double sensitivity;

        public double Sensitivity
        {
            get { return sensitivity; }
            set { sensitivity = value; }
        }


    }

    class RecordStruct
    {
        int seqcol;

        public int Seqcol
        {
            get { return seqcol; }
            set { seqcol = value; }
        }
        int conscol;

        public int Conscol
        {
            get { return conscol; }
            set { conscol = value; }
        }
        int poscol;

        public int Poscol
        {
            get { return poscol; }
            set { poscol = value; }
        }
        char rcflag;

        public char Rcflag
        {
            get { return rcflag; }
            set { rcflag = value; }
        }

    }
}
