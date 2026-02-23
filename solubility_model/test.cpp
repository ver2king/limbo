#include "opt_func/opt_func.hpp"

using namespace DATA;
using namespace MODEL;
using namespace OPT_FUNCTION;

namespace Test
{   
    using namespace DATA;
    using namespace MODEL;
    using namespace OPT_FUNCTION;

    //
    std::ifstream dataFile("/home/v183p176/Desktop/limbo/solubility_model/xH2_1.txt");
    PROPERTY dataProp = PROPERTY::H2LIQUID;
    String dataNamePrefix = "h2_liquid";
    //
    DATA::ExperimentalSolubilityData expSolData( dataFile, dataProp );
};

int main()
{   
    using namespace Test;
    expSolData.setPropertyNamePrefix( dataNamePrefix );
    expSolData.getPropertyData();
    //
    DataIdentifier sampleDataIdentifier = { "h2_liquid_4", dataProp };
    std::map < DataIdentifier, Data > expSolDataMap = expSolData._experimentalSolubilityData;
    Data expSolDataObj =  expSolDataMap[ sampleDataIdentifier ];
    //
    PrintTensor1D<double>( expSolDataObj._pressureData );
    PrintTensor1D<double>( expSolDataObj._temperatureData);
    PrintTensor2D<double>( expSolDataObj._propertyData );
    return 0;
};