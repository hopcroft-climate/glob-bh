class GST
    { 
      public:
          int set_history(vector<double> GST, int d);
          vector<double> get_history(); 
          int get_size();
          
      private:
		  vector<double> GST_history;
		           
    };

//================================================================================================    
int GST::set_history(vector<double> GST, int d)
{
	
	GST_history.resize(d);
	GST_history = GST;
	
	return 0;
}
//================================================================================================
vector<double> GST::get_history()
{
	vector<double>::iterator q;
	//for (q=GST_history.begin(); q!=GST_history.end(); q++)
	//{
	//	cerr << *q << "\t";
	//}
	//cerr << endl;
		
	return GST_history;
}
//================================================================================================
int GST::get_size()
{
	int i=0;
	i = GST_history.size();
	//cerr << "GST length" << i << endl;
	return i;
}
//================================================================================================
