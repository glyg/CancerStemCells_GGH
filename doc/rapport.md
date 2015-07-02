# Mise en place d'une simulation GGH de la prollifération dans les Neuroblastes

## Position du problème


Dans l'hypothèse des cellules souches cancéreuses (Cancer Stem Cells - CSC), la
population tumorale n'est pas homogène. Seule une sous-population, les CSC,
maintiennent la croissance de tumeur du fait de leur fort pouvoir prolifératif.
Cependant, une fraction des CSC peuvent se différencier en cellules cancéreuses
avec un pouvoir prolifératif plus faible.

Dans les tumeurs observées dans le modèle dévelopé par C. Maurange, les CSC
s'organisent en agrégats entourés de cellules différenciées. Pour saisir
finement l'origine et l'importance dans la prolifération de cette organisation
spatiale, un modèle biophysique peut être développé.

C. Maurange a fait appel à la société DamCB pour la mise en place d'une
simulation de la prolifération tumorale dans les neuroblastes de Drosophile.

Le travail c'est déroulé en 3 étapes:

* Identification du type de modélisation adaptée
* Installation/choix d'un logiciel de modélisation
* Première simulation à partir d'hypothèses simples

## Solution technique


Les modèles de Potts cellulaires, ou modèles de Glaziers Gardner Hoggs (GGH)
sont bien adaptés aux problématiques de prolifération et de différentiation
cellulaire. En effet, ils permettent de bien rendre compte des intéractions
cellule-cellule, qu'elles soient **biomécaniques** (adhésion) ou
**biochimiques** (signaux).


## Logiciel de simulation GGH: CompuCell3D

CompuCell3D has been developed exactly with the GGH framework, by the framework
authors themselves.

### Current Model specification in CompuCell3D:

Le modèle défini trois types cellulaires:

* The `CancerStemCell` type, representing the proliferative cell line
* The `Differentiated` type, less proliferative
* The `Medium` itself is defined as a cell type

#### Common properties

Those properties are stored in `../Sim2/Simulation/sim2.xml`


#### Xml specification


```xml
  <!-- Containr for the hole specificatiion of CPM (GGH) algorithm -->
     <Plugin Name="Contact">
        <!-- Specification of adhesion energies -->
        <Energy Type1="Medium" Type2="Medium">10.0</Energy>
        <Energy Type1="Medium" Type2="CancerStemCell">10.0</Energy>
        <Energy Type1="Medium" Type2="Differentiated">10.0</Energy>
        <Energy Type1="CancerStemCell" Type2="CancerStemCell">10.0</Energy>
        <Energy Type1="CancerStemCell" Type2="Differentiated">10.0</Energy>
        <Energy Type1="Differentiated" Type2="Differentiated">10.0</Energy>
        <NeighborOrder>1</NeighborOrder>
      </Plugin>
```


## Premiers travaux de modélisation


## Code Python

## Suites possibles



[1] Swat, M.H. et al. _Multi-Scale Modeling of Tissues Using CompuCell3D_.
Method Cell Biol 110, 325-366 (2012).
