# Mise en place d'une simulation GGH de la prollifération dans les Neuroblastes

## Position du problème


Dans l'hypothèse des cellules souches cancéreuses (Cancer Stem Cells - CSC), la
population tumorale n'est pas homogène. Seule une sous-population, les CSC,
maintiennent la croissance de tumeur du fait de leur fort pouvoir prolifératif.
Cependant, une fraction des CSC peuvent se différencier en cellules non cancéreuses (Non Cancerous Progenitor - NCP) avec un pouvoir prolifératif plus faible.

Dans les tumeurs observées dans le modèle dévelopé par C. Maurange, les CSC
s'organisent en agrégats entourés de NCP. Pour saisir
finnement l'origine et l'importance dans la prolifération de cette organisation
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

Le modèle GGH est décrit sur une grille de pixels fixe (Fig. 1 A). À chaque
pixel est associé un type (ci-dessous: vert, blanc ou bleu). À chaque pas de la
simulation (correspondant à un pas de temps discret), chaque pixel du modèle
peut changer d'état, selon un processus stochastique dépendant de l'énergie
associée aux intéractions avec ses voisins. Les changements de type favorables
énergétiquement sont privilégiés (Fig. 1 B).




![Schéma de principe du modèle de Potts cellulaire employé. ](../images/figure1.png)



### Logiciel de simulation GGH: CompuCell3D

CompuCell3D has been developed exactly with the GGH framework, by the framework
authors themselves.

### Current Model specification in CompuCell3D:

The model defines three possible pixel `Types`:

1. The `CancerStemCell` type, representing the proliferative cell line (CSC)
2. The `NonCancerous` type, less proliferative (NCP)
3. The surrounding `Medium`




#### Common properties

Those properties are stored in `../Sim2/Simulation/sim2.xml`


#### Xml specification


```xml
  <!-- Containr for the whole specificatiion of CPM (GGH) algorithm -->
     <Plugin Name="Contact">
        <!-- Specification of adhesion energies -->
        <Energy Type1="Medium" Type2="Medium">10.0</Energy>
        <Energy Type1="Medium" Type2="CancerStemCell">10.0</Energy>
        <Energy Type1="Medium" Type2="NonCancerous">10.0</Energy>
        <Energy Type1="CancerStemCell" Type2="CancerStemCell">10.0</Energy>
        <Energy Type1="CancerStemCell" Type2="NonCancerous">10.0</Energy>
        <Energy Type1="NonCancerous" Type2="NonCancerous">10.0</Energy>
        <NeighborOrder>1</NeighborOrder>
      </Plugin>
```


## Premiers travaux de modélisation


## Code Python

## Suites possibles



[1] Swat, M.H. et al. _Multi-Scale Modeling of Tissues Using CompuCell3D_.
Method Cell Biol 110, 325-366 (2012).
