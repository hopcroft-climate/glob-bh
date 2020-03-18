// labelling class for Voronoi tessellations taking account of degeneracy
// 10.7.07

class Label
    { 
      public:
          int set_label(vector<int> GST_id, int d);
          vector<int> get_label(); 
          int get_size();
          vector<int> get_code(int d);
          
      private:
		  vector<int> GST_label;
		           
    };

//================================================================================================    
int Label::set_label(vector<int> GST_id, int d)
{
	GST_label.resize(d);
	GST_label = GST_id;
	
	return 0;
}
//================================================================================================
	
vector<int> Label::get_label()
{
	//vector<int>::iterator q;
	//for (q=GST_label.begin(); q!=GST_label.end(); q++)
	//{
	//	cerr << *q << "\t";
	//}
	//cerr << endl;
		
	return GST_label;
}
//================================================================================================
int Label::get_size()
{
	int i=0;
	i = GST_label.size();
	//cerr << "Label length " << i << endl;
	return i;
}
//================================================================================================
vector<int> Label::get_code(int d)
{
	//vector<int>::iterator r;
	vector<int> letters;
	vector<int> code;
	code.resize(d,0);
	bool En = false;
	int label = 1;
   // int r = 0;
    int RL= 0;
	letters.push_back(GST_label[0]);
	code[0] = 1;
	//cerr <<"letters size " << letters.size() << endl;
	for (int q=1; q<d;  q++)	
	{
		
		En =false;
		
		for (int r=0; r<(int(letters.size())); r++)
		{
			if (GST_label[q] == letters[r]) {En=true;}
		}
		if (En  ==true)
		// if we've had this letter already
		{	
			for (int i=0; i<q; i++)
			{
				if (GST_label[q] == GST_label[i])
				{     
					 RL = i;
				}
			}
			
			//cerr <<" r already " << r << endl;
			code[q] = code[RL]  ; //write the repeated label !!!!
			
		}
		else 
		// if this is a new letter label
		{	
			label ++;	// move on to next label
			code[q] =label;	// add this to the code
			letters.push_back(GST_label[q]);	 //add it to letters used.
		}
		
	}
	//cerr << "size code " << code.size() << endl;
	//cerr <<"size GST_label " << GST_label.size() << endl;
	return code;
}

