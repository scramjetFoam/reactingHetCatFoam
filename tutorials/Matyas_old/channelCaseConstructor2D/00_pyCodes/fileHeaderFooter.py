def retFileHeaderStr():
    """Return list of file header strings for blockMeshDict."""
    retStr = []
    retStr.append('/*--------------------------------*- C++ -*----------------------------------*\ \n')
    retStr.append('  =========                 |                                                   \n')
    retStr.append('  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox             \n')
    retStr.append('   \\\\    /   O peration     | Website:  https://openfoam.org                    \n')
    retStr.append('    \\\\  /    A nd           | Version:  10                                      \n')
    retStr.append('     \\\\/     M anipulation  |                                                   \n')
    retStr.append('\*---------------------------------------------------------------------------*/ \n')
    
    retStr.append('FoamFile \n')
    retStr.append('{ \n \t version \t 1.0; \n \t format \t ascii; \n')
    retStr.append(' \t class \t\t dictionary; \n \t object \t blockMeshDict; \n} \n')
    retStr.append('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n\n')
    
    return retStr

def retFileFooterStr():
    """Returns list of file footer strings."""
    retStr = []
    retStr.append('// ************************************************************************* //\n\n')
    return retStr