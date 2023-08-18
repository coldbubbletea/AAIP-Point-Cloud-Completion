#include<bits/stdc++.h>
using namespace std;

int main()
{
    ifstream real_emd_file("8192_REAL_EMD.txt");
    ifstream comparison_emd_file("emd_value_comparison.txt");
    ofstream MSE_comparison("MSE_comparison(magnify emd a thousand times).txt");
    double MSE_AAIP=0.0;
    double MSE_emd2=0.0;
    string indent;
    getline(comparison_emd_file,indent);
    getline(comparison_emd_file,indent);
    for(auto i=0;i<80;i++)
    {
        
        string real_id,comparison_id;
        double real_emd;
        double AAIP_emd;
        double emd2;
        real_emd_file>>real_id>>real_emd;
        comparison_emd_file>>comparison_id>>AAIP_emd>>emd2;
        real_emd*=1000;
        AAIP_emd*=1000;
        emd2*=1000;
        MSE_AAIP+=(real_emd-AAIP_emd)*(real_emd-AAIP_emd);
        MSE_emd2+=(real_emd-emd2)*(real_emd-emd2);
    }
    MSE_comparison<<"AAIP_MSE: "<<MSE_AAIP/80<<endl;
    MSE_comparison<<"emd2_MSE: "<<MSE_emd2/80<<endl;
    return 0;
}