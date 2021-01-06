
/// Ancestor class of every material Hamiltonian generator 
/**
 * Note that virtual function has litte cost ~ switch
 * @author Dasol Kim
 */
#pragma once

#include <string>
#include <array>
#include <tuple>
#include <libconfig.h++>
#include "constant.h"

enum class InitialValueType
{
    UniformValence,
    FermiDirac,
    Custom
};

class BaseMaterial
{
public:
    /// Number of bands
    int Nband;
    /// Number of valence bands - only used for indicating which band is lowest conduction band
    int Nval;

    /// virtual desctructor is used to prevent memory leak
    virtual ~BaseMaterial() {};
    
    /// Generate custom inital value.
    /**
     * Generate custom initial value of density matrix. Should be used only when you need strange initial value.
     * This is not pure virtual function. One can ignore this when implementing Material class.
     * @param[out] _dmstore Initial value of density matrix
     * @param[in] _kpoint kpoints = {kx, ky, kz}
     */
    virtual void GenInitialValue(complex *_dmstore, std::array<double, Ndim> _kpoint)
    {
        cerr << "Error: Try to using custom inital value without implementing GenInitialValue.";
        exit(EXIT_FAILURE);
    };

    /// Generate Hamiltonian matrix - tight binding matrix
    /**
     * @param[out] _hstore Hamiltonian matrix
     * @param[in] _kpoint kpoints = {kx, ky, kz}
     */
    virtual void GenHamiltonian(complex *_hstore, std::array<double, Ndim> _kpoint) = 0;

    /// Generate dipole matrix in Wannier basis
    /**
     * @param[out] _dstore dipole matrix
     * @param[in] _kpoint kpoints = {kx, ky, kz}
     */
    virtual void GenDipole(complex **_dstore, std::array<double, Ndim> _kpoint) = 0;

    /// Generate interband and intraband current matrix
    /**
     * Remember 
     * @param[out] _jstore store \f$ - \frac{\partial H(\textbf{k})}{\partial \textbf{k}} + i [\textbf{D}(\textbf{k}), H(\textbf{k})] \f$
     * @param[in] _kpoint kpoints = {kx, ky, kz}
     */
    virtual void GenJMatrix(complex **_jstore, std::array<double, Ndim> _kpoint) = 0;

    /// Generate Brillouin zone information
    /**
     * To use this function, see below example. @code{.cpp} auto [axes, origin] = GenBrillouinzone(); @endcode
     * @return Brillouin zone axes and origin points are returned as a tuple
     */
    virtual std::tuple<std::array<double, Ndim*Ndim>, std::array<double, Ndim> > GenBrillouinzone() = 0;

    virtual void PrintMaterialInformation() = 0;
};