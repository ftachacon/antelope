/// collection of material headers
/**
 * @author Dasol Kim
*/
#pragma once

#include "Materials/Haldane.h"
#include "Materials/Haldane2L.h"
#include "Materials/WilsonMass.h"
#include "Materials/Bi2Se3surf.h"
#include "Materials/TMDCs.h"
#include "Materials/KaneMele.h"
#include "Materials/XYHlattice.h"

//#include "Materials/Wannier90.h"


// global function in header file should be avoided, but we will use single cpp file anyway.
BaseMaterial* InitializeMaterial(std::string targetMaterial, const libconfig::Setting* _cfg)
{
    BaseMaterial *material;
    const libconfig::Setting &cfg = (*_cfg);
    if (targetMaterial == "Haldane")
    {
        material = new Haldane( &(cfg[targetMaterial.c_str()]) );
    }
    else if (targetMaterial == "Haldane2L")
    {
        material = new Haldane2L( &(cfg[targetMaterial.c_str()]) );
    }
    else if (targetMaterial == "KaneMele")
    {
        material = new KaneMele( &(cfg[targetMaterial.c_str()]) );
    }
    else if (targetMaterial == "WilsonMass")
    {
        material = new WilsonMass( &(cfg[targetMaterial.c_str()]) );
    }
    else if (targetMaterial == "Bi2Se3surf")
    {
        material = new BieSe3surf( &cfg );
    }
    else if (targetMaterial == "TMDCs")
    {
        material = new TMDCs( &cfg[targetMaterial.c_str()] );
    }
    else if (targetMaterial == "XYHlattice")
    {
        material = new XYHlattice( &cfg[targetMaterial.c_str() ]);
    }
    else
    {
        std::cerr << "Undefined Material\n";
        std::exit(EXIT_FAILURE);
        //material = new Wannier90( &(cfg[targetMaterial.c_str()]) );
        //isDipoleZero = dynamic_cast<Wannier90*>(material)->isDipoleZero;
        //vec_lattice = dynamic_cast<Wannier90*>(material)->vec_lattice;
        //isWannier90 = true;
    }
    return material;
}