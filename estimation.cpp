#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/helpers/ShortcutsGeometry.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/boards/Board3DTo2D.h"

#include <DGtal/topology/SCellsFunctors.h>
#include <DGtal/images/SimpleThresholdForegroundPredicate.h>

#include "DGtal/geometry/volumes/TangencyComputer.h"
#include "DGtal/geometry/surfaces/COBAGenericNaivePlaneComputer.h"
#include "DGtal/geometry/surfaces/ChordGenericNaivePlaneComputer.h"

#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"

#include "DGtal/geometry/surfaces/DigitalSurfacePredicate.h"
#include "DGtal/geometry/surfaces/estimation/PlaneProbingParallelepipedEstimator.h"
#include "DGtal/geometry/surfaces/estimation/PlaneProbingDigitalSurfaceLocalEstimator.h"

using namespace std;
using namespace DGtal;
using namespace Z3i;

// auteur : Aude Marêché

// usage :
// ./reset vol_file.vol abscisse ordonnée cote rayon


// intersection entre un ensemble set de points et le cube de côté d autour d'un point p
DigitalSet tricube(Point p, DigitalSet set, int d) {
    DigitalSet ret(set.domain());
    for (int x = p[0] - d; x < p[0] + d + 1; x++) {
        for (int y = p[1] - d; y < p[1] + d + 1; y++) {
            for (int z = p[2] - d; z < p[2] + d + 1; z++) {
                Point p = Point(x, y, z);
                if (set(p)) {
                    ret.insertNew(p);
    }   }   }   }
    return ret;
}


// fonction qui calcule le Trimmed Grassmann Average (moyenne) d'un vector de normales
// le résultat est normalisé
RealPoint TGA(vector<RealPoint> normales, float trim_per) {

    // initialisation de la moyenne selon une distribution uniforme
    static default_random_engine gen;
    static uniform_real_distribution<double> dist(-1.0,1.0);
    RealPoint q = RealPoint(dist(gen), dist(gen), dist(gen));

    // normalisation des vecteurs normaux aux secteurs
    for (int i = 0; i < normales.size(); i++) {
        normales[i] = normales[i].getNormalized();
    }

    // calcul de la moyenne en utilisant le Trimmed Grassmann Average (TGA)

    for (int i = 0; i < 10; i++) { // nb max d'itérations purement arbitraire (mais a priori suffisant dans nos cas) (voir la vraie condition d'arrêt (plus d'évolution de la moyenne) en fin de boucle si besoin)
        vector<RealPoint> signed_normales(normales);
        vector<double> dim1, dim2, dim3;

        // récupération des valeurs
        for (int n = 0; n < normales.size(); n++) {

            // isolement des différentes dimensions pour pouvoir faire la moyenne "pixel-wise" (-> robuste)
            dim1.push_back(signed_normales[n][0]);
            dim2.push_back(signed_normales[n][1]);
            dim3.push_back(signed_normales[n][2]);
        }


        // trimming
        // calcul des indices de la zone à garder
        size_t ind_deb = dim1.size() * (trim_per/100);
        size_t ind_fin = dim1.size() - ind_deb;

        // on s'assure qu'on garde au moins une valeur
        if (ind_fin - ind_deb < 1) {
            ind_fin += 1;
        }

        // tri pour pouvoir garder les énièmes éléments
        sort(dim1.begin(), dim1.end());
        sort(dim2.begin(), dim2.end());
        sort(dim3.begin(), dim3.end());

        // calcul de la moyenne courante
        int nb_normales = 0;
        RealPoint somme_normales = RealPoint(0, 0, 0);
        double somme_dim1 = 0;
        double somme_dim2 = 0;
        double somme_dim3 = 0;

        for (int ind = ind_deb; ind < ind_fin; ind++) {
            nb_normales += 1;

            somme_dim1 += dim1[ind];
            somme_dim2 += dim2[ind];
            somme_dim3 += dim3[ind];
        }

        RealPoint q_j = RealPoint(somme_dim1/nb_normales, somme_dim2/nb_normales, somme_dim3/nb_normales).getNormalized();

        // condition d'arrêt
        if (not (q_j - q).norm() > 0) {
            i = 20; // valeur arbitraire pour sortir de la boucle
        }
        q = q_j;
    }

    return q;
}



// fonction qui calcule le barycentre géodésique (moyenne de Fréchet) d'un vector de normales
// le résultat est normalisé
RealPoint BG(vector<RealPoint> normales) {

    // normalisation
    for (int i = 0; i < normales.size(); i++) {
        normales[i] = normales[i].getNormalized();
    }

    // version grid search/estimation

    // initialisation
    double pas = 0.01; // descends pas à 0.001, ça prend beaucoup trop de temps
    RealPoint moyenne = RealPoint(0, 0, 0);
    double phi_min = std::numeric_limits<double>::max();

    double theta_min = 0;
    double theta_max = M_PI;
    double delta_min = - M_PI;
    double delta_max = M_PI;


    // parcours de la grille
    for (double theta = min(theta_min, theta_max); theta <= max(theta_min, theta_max); theta += pas) {

        double pas_delta = pas; 
        if (theta != 0) {
            pas_delta = pas / sin(theta);
        }
        for (double delta = min(delta_min, delta_max); delta <= max(delta_min, delta_max); delta += pas_delta) {

            RealPoint p = RealPoint(sin(theta) * cos(delta), sin(theta) * sin(delta), cos(theta)).getNormalized();
            double phi = 0;

            for (int i = 0; i < normales.size(); i++) {
                double dot_prod = p.dot(normales[i]);
                if (dot_prod < 1) {
                    phi += acos(dot_prod) * acos(dot_prod);
                }
            }

            if (phi < phi_min) {
                phi_min = phi;
                moyenne = p;
            }
        }
    }
    return moyenne.getNormalized();
}




// main -------------------------------------------------------------------
int main(int argc, char** argv) {

    typedef Shortcuts<Z3i::KSpace> SH3;
    typedef ShortcutsGeometry<Z3i::KSpace> SHG3;
    auto params = SH3::defaultParameters();
    typedef KSpace::SurfelSet SurfelSet;

    // chargement de l'objet (.vol)
    cout << "\nInitialisation...\n";
    auto binary_image = SH3::makeBinaryImage(argv[1], params);
    auto domain = (*(binary_image.get())).domain();
    auto K = SH3::getKSpace(binary_image);
    auto surface = SH3::makeDigitalSurface(binary_image, K, params);
    auto surfels = SH3::getSurfelRange(surface, params);
    cout << "Objet : " << argv[1] << "\n";

    // récupération du point considéré et du rayon maximal
    Point p(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    cout << "Point considéré : " << p << "\n\n";
    int rayon = atoi(argv[5]);

    // map surfel -> voxel
    functors::SCellToInnerPoint<KSpace> map(K);


    // récupération des surfels de la zone de travail (r + 2)
    L2Metric norme_l2;
    SurfelSet surfels_intersection;
    cout << "Récupération des surfels de la zone autour du point " << p << "...\n";
    for (SCell s : surfels) {
        if (norme_l2(p, map(s)) <= rayon + 2) {
            surfels_intersection.insert(s);
        }
    }

    // calcul des "points de la bordure intérieure" pour le Tangency Computer
    vector<Point> points_int;
    std::map<Point, int> point_to_id;
    int id = 0;

    // parcours des surfels
    cout << "Calcul des points internes...\n";
    for (auto s : surfels_intersection) {
        Point q = map(s);
        // si on n'a pas encore vu le point
        if (point_to_id.find(q) == point_to_id.end()) {
            points_int.push_back(q);
            point_to_id[q] = id ++;
        }
    }
    // récupération de l'indice de p
    int id_p = point_to_id[p];


    // calcul de la carte des distances
    cout << "Calcul de la carte des distances...\n\n";
    typedef TangencyComputer<KSpace>::Index Index;
    TangencyComputer<KSpace> TC(K);
    TC.init(points_int.cbegin(), points_int.cend());
    auto SP = TC.makeShortestPaths(sqrt(3.0));
    SP.init(id_p);
    double last_distance = 0.0;

    while (!SP.finished()) {
        last_distance = std::get<2>(SP.current());
        SP.expand();
    }


    // recherche du rayon minimum où on n'est plus planaire
    COBAGenericNaivePlaneComputer<Z3, int64_t> plan_max;
    plan_max.init(100, 1, 1);
    DigitalSet voxels_intersection(domain);
    int petit_rayon = 1;
    bool planaire = true;

    while (petit_rayon <= rayon && planaire) {
        for (SCell s : surfels_intersection) {
            if (SP.distance(point_to_id[map(s)]) <= petit_rayon) {
                voxels_intersection.insertNew(map(s));
            }
        }
        if (not plan_max.extend(voxels_intersection.begin(), voxels_intersection.end())) {
            planaire = false;
        }
        petit_rayon ++;
    }
    rayon = petit_rayon - 1;
    cout << "Premier rayon d'intersection non planaire : " << rayon << "\n\n";


    // définition de la color map par rapport au rayon
    GradientColorMap<long> gradient(0, rayon);
    gradient.addColor(Color(255, 0, 128));
    gradient.addColor(Color::Magenta);
    gradient.addColor(Color::Cyan);
    gradient.addColor(Color::Yellow);
    gradient.addColor(Color(255, 128, 0));
    gradient.addColor(Color::Red);
    gradient.addColor(Color::Black);


    // calcul des points de la bordure de l'intersection (i.e. les points dont au moins un surfel a un voisin extérieur à l'intersection en 4-connexité)
    cout << "Calcul de la bordure/contour de l'intersection...\n";
    DigitalSet bordure(domain);
    std::map<Cell, set<Cell>> pointels; // lignels to pointels
    std::map<Cell, set<Cell>> lignels; // pointels to lignels
    std::map<Cell, Point> voxel; // lignels to voxel
    Cell dep = K.uCell(Point(0, 0, 0));

    for (SCell s : surfels_intersection) {
        if (SP.distance(point_to_id[map(s)]) <= rayon) {

            // récupération des voisins de s
            SH3::SurfelRange vois;
            back_insert_iterator<SH3::SurfelRange> writeIt = back_inserter(vois);
            (*(surface.get())).writeNeighbors(writeIt, s);

            // on regarde si on a un voisin en dehors de l'intersection (on est au bord)
            bool is_bord = false;
            for (SCell v : vois) {
                if (SP.distance(point_to_id[map(v)]) > rayon) {
                    is_bord = true;
                }
            }
            if (is_bord) {
                bordure.insertNew(map(s));

                // si on est au bord, on récupère les voxels, lignels et pointels dont on a besoin pour trier les points dans l'ordre
                Cells lignels_s = K.uLowerIncident(K.unsigns(s));
                for (Cell ls : lignels_s) {
                    voxel[ls] = map(s);
                }

                for (SCell v : vois) {
                    if (SP.distance(point_to_id[map(v)]) > rayon) {

                        Cells lignels_v = K.uLowerIncident(K.unsigns(v));
                        for (Cell l : lignels_v) {
                            for (Cell ls : lignels_s) {

                                // le lignel est au bord
                                if (ls == l) {
                                    // point de départ de la bordure triée
                                    if (dep == K.uCell(Point(0, 0, 0))) {
                                        dep = l;
                                    }
                                    Cells pointels_l =  K.uLowerIncident(l);
                                    for (Cell p : pointels_l) {
                                        pointels[l].insert(p);
                                        lignels[p].insert(l);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "Nombre de point sur la bordure : " << bordure.size() << "\n";


    // tri des points de la bordure
    vector<Cell> lignels_bordure;
    lignels_bordure.push_back(dep);

    auto it = pointels[dep].begin();
    Point Kp = K.uKCoords(K.uSpel(p));

    // détermination de la direction de parcours et du premier voisin

    // récupération de la normale à la surface au point p
    Vector N = Vector(0, 0, 0);
    for (SCell s : surfels_intersection) {
        if (map(s) == p) {
            auto j = K.sOrthDir(s);
            bool d = K.sDirect(s, j);
            N += K.sCoords(K.sIncident(s, j, !d)) - K.sCoords(K.sIncident(s, j, d));
        }
    }

    Point B = K.uKCoords(dep) - Kp;
    Point T1 = K.uKCoords(*it) - Kp;
    Cell cur_pointel = *it;
    Point T2 = K.uKCoords(*(++it)) - Kp;

    // choix du pointel dans la bonne direction
    Point prod = B.crossProduct(T1);
    if (prod.dot(N) < 0) {
        prod = B.crossProduct(T2);
        cur_pointel = *it;
    }
    if (prod.dot(N) < 0) {
        cout << "oups ?\n";
        exit(1);
    }

    // premier voisin
    pointels[dep].erase(cur_pointel);
    it = lignels[cur_pointel].begin();
    Cell cur_lignel = *it;
    if (cur_lignel == dep) {
        cur_lignel = *(++it);
    }
    lignels[cur_pointel].erase(cur_lignel);
    pointels[cur_lignel].erase(cur_pointel);
    cur_pointel = *(pointels[cur_lignel].begin());

    // parcours des lignels et pointels du bord
    int cpt = 0;
    while (cur_lignel != dep && cpt < 400) {
        lignels_bordure.push_back(cur_lignel);
        lignels[cur_pointel].erase(cur_lignel);
        cur_lignel = *(lignels[cur_pointel].begin());
        pointels[cur_lignel].erase(cur_pointel);
        cur_pointel = *(pointels[cur_lignel].begin());
        cpt ++;
    }
    if (cpt == 400) {
        cout << "on ne peut pas compléter le contour :-/\n"; // on n'a pas réussi à parcourir toute la bordure correctement
    }

    // récupération des voxels dans l'ordre
    vector<Point> bord_triee;
    Point cur_voxel = voxel[dep];
    for (Cell l : lignels_bordure) {
        if (voxel[l] != cur_voxel) {
            cur_voxel = voxel[l];
            bord_triee.push_back(cur_voxel);
        }
    }
    if (bord_triee[0] != voxel[dep] && cur_voxel != voxel[dep]) {
        bord_triee.push_back(voxel[dep]);
    }

    cout << "Nombre de points triés : " << bord_triee.size() << "\n\n";


    // calcul des plus courts chemins (tangents) (paths) entre les points de la bordure et el centre p
    cout << "Calcul des plus courts chemins tangents entre les points de la bordure et le centre " << p << "...\n";
    vector<Index> dest = {id_p};
    vector<Index> sources;
    for (Point q : bordure) {
        sources.push_back(point_to_id[q]);
    }
    auto paths = TC.shortestPaths(sources, dest, sqrt(3.0));


    // parcours des chemins
    cout << "Calcul des chemins en voxels...\n\n";
    vector<vector<Point>> chemins; // si ça marche on enlèvera le "chemin" intermédiaire -> ou pas
    COBAGenericNaivePlaneComputer<Z3, int64_t> plan;
    plan.init(100, 1, 1);
    vector<int> plans(5, 0);
    vector<int> pasplans(5, 0);

    // parcours des points de la bordure
    for (int nb = 0; nb < paths.size(); nb ++) {
        vector<Point> chemin;

        // calcul des chemins de voxels à partir des paths (segments de droites réelles/lignes brisées/polylignes)
        for (int i = 0; i < paths[nb].size() - 1; i ++) {

            // affichage du segment de droite réelle
            Point p1 = TC.point(paths[nb][i]);
            Point p2 = TC.point(paths[nb][i + 1]);
            Point droite = p2 - p1;

            // initialisation du parcours
            Point cur_voxel = p1;
            Point prev_voxel = cur_voxel;
            Point prev_prev_voxel = cur_voxel;
            chemin.push_back(cur_voxel);


            // parcours des voxels
            while (cur_voxel != p2) {

                // récupération des voisins
                double angle_min = 7;
                Point voxel_min = cur_voxel;
                DigitalSet vois = tricube(cur_voxel, voxels_intersection, 1);
                vois.erase(cur_voxel);
                vois.erase(prev_voxel);
                vois.erase(prev_prev_voxel);

                // calcul du voisin le plus proche (par l'angle) de la droite réelle
                for (Point v : vois) {
                    Point candidat = v - p1;
                    double angle;
                    double produit = droite.dot(candidat) / (droite.norm() * candidat.norm());
                    produit >= 1 ? angle = 0 : angle = acos(produit);

                    if (angle < angle_min && SP.distance(point_to_id[v]) <= SP.distance(point_to_id[cur_voxel])) {
                        angle_min = angle;
                        voxel_min = v;
                    }
                }

                // on passe au point suivant
                prev_prev_voxel = prev_voxel;
                prev_voxel = cur_voxel;
                cur_voxel = voxel_min;
                chemin.push_back(cur_voxel);
            }

            // calcul de l'angle de la ligne brisée
            if (paths[nb].size() > 2 && i > 0) {
                Point p0 = TC.point(paths[nb][i - 1]);
                Vector u = p0 - p1;
                Vector v = p2 - p1;
            }

        }

        // test de planéité
        if (plan.isExtendable(chemin.begin(), chemin.end())) {
            plans[paths[nb].size() - 2] += 1;
            chemins.push_back(chemin);
        } else {
            pasplans[paths[nb].size() - 2] += 1;
        }
    }

    // affichage des statistiques sur la planéité des chemins calculés
    cout << "Nombre de brisures :\t    0  1 2 3 4\n";
    cout << "Plans pour x brisures :\t    ";
    for (int i = 0 ; i < plans.size(); i ++) {
        cout << plans[i] << " ";
    }
    cout << "\nNon-plans pour x brisures : ";
    for (int i = 0 ; i < pasplans.size(); i ++) {
        cout << pasplans[i] << " ";
    }
    cout << "\n\n\n";


    // tri des chemins dans l'ordre de la bordure
    vector<DigitalSet> chemins_tries;
    for (int i = 0; i < bord_triee.size(); i++) {
        int j = 0;
        while (j < chemins.size() && chemins[j][0] != bord_triee[i]) {
            j ++;
        }
        if (j < chemins.size()) {
            DigitalSet chemin(domain);
            chemin.insert(chemins[j].begin(), chemins[j].end());
            chemins_tries.push_back(chemin);
            chemins.erase(chemins.begin() + j);
        }
    }

    // récupérations des points n'appartenant pas à un chemin
    DigitalSet points_seuls(domain);
    for (Point v : voxels_intersection) {
        bool in = false;
        int i = 0;
        while (not in && i < chemins_tries.size()) {
            in = chemins_tries[i](v);
            i ++;
        }
        if (not in) {
            points_seuls.insert(v);
        }
    }


    // construction des secteurs et estimation des normales
    cout << "Construction des secteurs planaires...\n";

    // affichage
    QApplication application(argc, argv);
    Viewer3D<> viewer;
    viewer.show();

    viewer.setFillTransparency(176);
    viewer << CustomColors3D(Color::Black, Color::Black);
    viewer << p;

    vector<Color> colors = {Color::Red, Color::Yellow, Color::Green, Color::Cyan, Color(0, 170, 255), Color::Blue, Color(170, 0, 255), Color::Magenta, Color(255, 0, 127), Color(170, 0, 0)};


    // calcul des différents découpages
    vector<RealPoint> normales_COBA_decoups_TGA;
    vector<RealPoint> normales_Chord_decoups_TGA;
    vector<RealPoint> normales_COBA_decoups_BG;
    vector<RealPoint> normales_Chord_decoups_BG;
    int nb_secteurs_min = bordure.size();
    int nb_secteurs_max = 0;
    int nb_secteurs_moy = 0;
    int nb_points_seuls_min = voxels_intersection.size();
    int nb_points_seuls_max = 0;
    int nb_points_seuls_moy = 0;

    // parcours des points de la bordure
    for (int dep = 0; dep < chemins_tries.size(); dep ++) {
        cout << "Découpage n°" << dep << "\n";
        vector<COBAGenericNaivePlaneComputer<Z3, int64_t>> secteurs_COBA;
        vector<ChordGenericNaivePlaneComputer<Space, Point, int64_t>> secteurs_Chord;
        vector<RealPoint> normales_COBA_secteurs;
        vector<RealPoint> normales_Chord_secteurs;
        DigitalSet p_s(points_seuls);

        COBAGenericNaivePlaneComputer<Z3, int64_t> secteur_COBA;
        secteur_COBA.init(100, 1, 1);
        ChordGenericNaivePlaneComputer<Space, Point, int64_t> secteur_Chord;
        secteur_Chord.init(1, 1);
    

        // ajout des chemins
        for (int i = dep; i < chemins_tries.size() + dep; i ++) {

            if (not secteur_COBA.extend(chemins_tries[i%chemins_tries.size()].begin(), chemins_tries[i%chemins_tries.size()].end())) { // si on n'est plus planaire
                DigitalSet ps(p_s);
                for (Point q : p_s) {
                    DigitalSet secteur(domain);
                    secteur.insert(secteur_COBA.begin(), secteur_COBA.end());
                    if (tricube(q, secteur, 1).size() > 2) {
                        if (secteur_COBA.extend(q)) {
                            ps.erase(q);
                        }
                    }
                }
                p_s = ps;

                secteurs_COBA.push_back(secteur_COBA);

                // estimations de la normale au secteur
                // COBA
                RealPoint normale;
                secteur_COBA.getNormal(normale);
                normales_COBA_secteurs.push_back(normale);
                //cout << "Normale estimée par COBA : " << normale << "\n";

                // Chord
                secteur_Chord.extend(secteur_COBA.begin(), secteur_COBA.end());
                secteurs_Chord.push_back(secteur_Chord);
                secteur_Chord.getNormal(normale);
                normale = normale.getNormalized();
                normales_Chord_secteurs.push_back(normale);
                //cout << "Normale estimée par Chord : " << normale << "\n";

                // ré-initialisation pour le prochain secteur
                secteur_COBA.init(100, 1, 1);
                secteur_Chord.init(1, 1);
                if (not secteur_COBA.extend(chemins_tries[i%chemins_tries.size()].begin(), chemins_tries[i%chemins_tries.size()].end())) {
                    cout << "Un chemin non planaire ? Bouh :-(\n"; // normalement il ne devrait pas y en avoir :-)
                    exit(1);
                }
            }  else { // on est encore planaire
                // on essaie d'ajouter les points isolés qui ont au moins trois voisins dans le secteur
                DigitalSet ps(p_s);
                for (Point q : p_s) {
                    DigitalSet secteur(domain);
                    secteur.insert(secteur_COBA.begin(), secteur_COBA.end());
                    if (tricube(q, secteur, 1).size() > 2) {
                        if (secteur_COBA.extend(q)) {
                            ps.erase(q);
                        }
                    }
                }
                p_s = ps;
            }
            
        }

        // on essaie d'ajouter les points encore isolés qui ont au moins trois voisins dans le dernier secteur
        DigitalSet ps(p_s);
        for (Point q : p_s) {
            DigitalSet secteur(domain);
            secteur.insert(secteur_COBA.begin(), secteur_COBA.end());
            if (tricube(q, secteur, 1).size() > 2) {
                if (secteur_COBA.extend(q)) {
                    ps.erase(q);
                }
            }
        }
        p_s = ps;

        // on essaie de fusionner le premier et le dernier secteur
        if (secteurs_COBA.size() > 0 && secteur_COBA.extend(secteurs_COBA[0].begin(), secteurs_COBA[0].end())) {
            secteurs_COBA[0] = secteur_COBA;

            // ré-estimation de la normale du premier secteur
            // COBA
            RealPoint normale;
            secteur_COBA.getNormal(normale);
            normales_COBA_secteurs[0] = normale;
            //cout << "Normale estimée par COBA : " << normale << "\n";

            // Chord
            secteur_Chord.extend(secteur_COBA.begin(), secteur_COBA.end());
            secteurs_Chord[0] = secteur_Chord;
            secteur_Chord.getNormal(normale);
            normale = normale.getNormalized();
            normales_Chord_secteurs[0] = normale;
            //cout << "Normale estimée par Chord : " << normale << "\n\n";

        } else {
            secteurs_COBA.push_back(secteur_COBA);

            // estimation de la normale au dernier secteur
            // COBA
            RealPoint normale;
            secteur_COBA.getNormal(normale);
            normales_COBA_secteurs.push_back(normale);
            //cout << "Normale estimée par COBA : " << normale << "\n";

            // Chord
            secteur_Chord.extend(secteur_COBA.begin(), secteur_COBA.end());
            secteurs_Chord.push_back(secteur_Chord);
            secteur_Chord.getNormal(normale);
            normale = normale.getNormalized();
            normales_Chord_secteurs.push_back(normale);
            //cout << "Normale estimée par Chord : " << normale << "\n\n";
        }


        // calcul des normales moyennes
        vector<RealPoint> normales_proport_COBA; // avec pondération par la taille des secteurs
        vector<RealPoint> normales_proport_Chord;
        vector<RealPoint> normales_pas_proport_COBA; // sans pondération par la taille des secteurs
        vector<RealPoint> normales_pas_proport_Chord;
        for (int i = 0; i < secteurs_COBA.size(); i ++) {
            RealPoint nor;
            secteurs_Chord[i].getNormal(nor);
            nor = nor.getNormalized();
            for (Point s : secteurs_COBA[i]) {
                normales_proport_COBA.push_back(normales_COBA_secteurs[i]);
                 normales_proport_Chord.push_back(normales_Chord_secteurs[i]);
            }
             normales_pas_proport_COBA.push_back(normales_COBA_secteurs[i]);
             normales_pas_proport_Chord.push_back(normales_Chord_secteurs[i]);
        }


        // affichage de "statistiques" sur les secteurs
        cout << "Nombre de secteurs : " << secteurs_COBA.size() << "\n";
        if (secteurs_COBA.size() < nb_secteurs_min) {
            nb_secteurs_min = secteurs_COBA.size();
        } else if (secteurs_COBA.size() > nb_secteurs_max) {
            nb_secteurs_max = secteurs_COBA.size();
        }
        nb_secteurs_moy += secteurs_COBA.size();

        cout << "Tailles des secteurs : ";
        for (auto s : secteurs_COBA) {
            cout << s.size() << " ";
        }
        cout << "; points hors secteurs : " << p_s.size() << "\n";
        if (p_s.size() < nb_points_seuls_min) {
            nb_points_seuls_min = p_s.size();
        } else if (p_s.size() > nb_points_seuls_max) {
            nb_points_seuls_max = p_s.size();
        }
        nb_points_seuls_moy += p_s.size();


        // TGA (Trimmed Grassmann Average)
        RealPoint moy = TGA(normales_proport_COBA, 50);
        normales_COBA_decoups_TGA.push_back(moy);
        cout << "\nnormale moyenne TGA 50 COBA : " << moy << "\n";
        moy = TGA(normales_proport_Chord, 50);
        normales_Chord_decoups_TGA.push_back(moy);
        cout << "normale moyenne TGA 50 Chord : " << moy << "\n";

        // barycentre géodésique
        moy = BG(normales_proport_COBA);
        normales_COBA_decoups_BG.push_back(moy);
        cout << "normale moyenne BG COBA : " << moy << "\n";
        moy = BG(normales_proport_Chord);
        normales_Chord_decoups_BG.push_back(moy);
        cout << "normale moyenne BG Chord : " << moy << "\n\n";


        // affichage du dernier découpage à titre d'exemple
        if (dep == chemins_tries.size() - 1) {
            // pour enregistrer une image
            /*
            Board3DTo2D<Space, KSpace> board;
            board << domain;
            board << CameraPosition(-14, 17, 60)
                  << CameraDirection(-0.32681, 0.216452, -0.919969)
                  << CameraUpVector(0, 1, 0);
            board << CameraZNearFar(0.57, 163);

            board << SetMode3D(p.className(), "Paving");
            board << CustomColors3D(Color::Black, Color::Black);
            board << p;
            */

            // ajout des points de l'intersection
            for (int i = 0; i < secteurs_COBA.size(); i++) {
                viewer << CustomColors3D(colors[i%colors.size()], colors[i%colors.size()]);
                //board << CustomColors3D(colors[i%colors.size()], colors[i%colors.size()]);
                for (Point q : secteurs_COBA[i]) {
                    viewer << q;
                    //board << q;
                }
            }

            // affichage des autres voxels de l'intersection
            
            for (auto s : surfels) {
                double dist = SP.distance(point_to_id[map(s)]);
                if (voxels_intersection(map(s))) {
                    // points hors secteurs
                    viewer << CustomColors3D(Color::Gray, Color::Gray);
                    //board << CustomColors3D(Color::Gray, Color::Gray);

                    // pour n'ajouter au board QUE les points vraiment hors secteurs
                    /*
                    bool out = false;
                    for (Point q : p_s) {
                        if (q == map(s)) {
                            out = true;
                        }
                    }
                    if (out) {
                        board << map(s);
                    }
                    */

                } else {
                    viewer << CustomColors3D(Color::White, Color::White);
                }
                viewer << map(s);
            }
            viewer << Viewer3D<>::updateDisplay;
            //application.exec();

            /*
            if (dep < 10) {
                board.saveCairo(("test/secteurs_0" + to_string(dep) + ".png").c_str(), Board3DTo2D<Space, KSpace>::CairoPNG, 600*2, 400*2);
            } else {
                board.saveCairo(("test/secteurs_" + to_string(dep) + ".png").c_str(), Board3DTo2D<Space, KSpace>::CairoPNG, 600*2, 400*2);
            }
            */
        }
    }
    cout << "\n";


    // statistiques sur le nombre de secteurs
    cout << "Plus petit nombre de secteurs : " << nb_secteurs_min << "\n";
    cout << "Plus grand nombre de secteurs : " << nb_secteurs_max << "\n";
    cout << "Nombre moyen de secteurs : " << double(nb_secteurs_moy) / chemins_tries.size() << "\n";
    cout << "Plus petit nombre de points isolés : " << nb_points_seuls_min << "\n";
    cout << "Plus grand nombre de points isolés : " << nb_points_seuls_max << "\n";
    cout << "Nombre moyen de points isolés : " << double(nb_points_seuls_moy) / chemins_tries.size() << "\n\n\n";


    // affichage sur la sphère
    Viewer3D<> viewer_sphere;
    viewer_sphere.show();
    viewer_sphere << CustomColors3D(Color::Gray, Color::Gray);
    viewer_sphere.addBall(RealPoint(0, 0, 0), 1, 100);

    //affichage de l'echelle
    viewer_sphere << CustomColors3D(Color::Black, Color::Black);
    for (int i = 0; i < 360; i++) {
        double rad = i * M_PI / 180;
        viewer_sphere.addBall(RealPoint(0, cos(rad), sin(rad)), 0.002);
        viewer_sphere.addBall(RealPoint(cos(rad), 0, sin(rad)), 0.002);
        viewer_sphere.addBall(RealPoint(cos(rad), sin(rad), 0), 0.002);
    }


    // calcul et affichage des moyennes sur les découpages et des moyennes globales
    cout << "Normales moyennes sur tous les découpages :\n";
    // TGA 50 COBA
    viewer_sphere << CustomColors3D(Color::Red, Color::Red);
    for (RealPoint n : normales_COBA_decoups_TGA) {
        viewer_sphere.addBall(n, 0.002);
    }
    RealPoint moy = TGA(normales_COBA_decoups_TGA, 50);
    cout << "TGA 50 COBA (rouge) : " << moy << "\n";
    viewer_sphere.addBall(moy, 0.004);

    // TGA 50 Chord
    viewer_sphere << CustomColors3D(Color::Green, Color::Green);
    for (RealPoint n : normales_Chord_decoups_TGA) {
        viewer_sphere.addBall(n, 0.002);
    }
    moy = TGA(normales_Chord_decoups_TGA, 50);
    cout << "TGA 50 Chord (vert) : " << moy << "\n";
    viewer_sphere.addBall(moy, 0.004);

    // BG COBA
    viewer_sphere << CustomColors3D(Color::Yellow, Color::Yellow);
    for (RealPoint n : normales_COBA_decoups_BG) {
        viewer_sphere.addBall(n, 0.002);
    }
    moy = BG(normales_COBA_decoups_BG);
    cout << "BG COBA (jaune) : " << moy << "\n";
    viewer_sphere.addBall(moy, 0.004);

    // BG Chord
    viewer_sphere << CustomColors3D(Color::Cyan, Color::Cyan);
    for (RealPoint n : normales_Chord_decoups_BG) {
        viewer_sphere.addBall(n, 0.002);
    }
    moy = BG(normales_Chord_decoups_BG);
    cout << "BG Chord (cyan) : " << moy << "\n\n";
    viewer_sphere.addBall(moy, 0.004);


    // estimation de la normale par Integral Invariant
    viewer_sphere << CustomColors3D(Color::Blue, Color::Blue);
    cout << "Estimation de la normale par Integral Invariant (bleu) :\n";

    double h = 0.25;
    typedef ImageSelector< Z3i::Domain, bool >::Type Image;
    typedef functors::SimpleThresholdForegroundPredicate< Image > ImagePredicate;
    Image image = VolReader<Image>::importVol(argv[1]);
    ImagePredicate predicate = ImagePredicate(image, 0);

    functors::IINormalDirectionFunctor<Space> normalFunctor;
    IntegralInvariantCovarianceEstimator<KSpace, ImagePredicate, functors::IINormalDirectionFunctor<Space>> normalEstimator(normalFunctor);
    normalEstimator.attach(K, predicate);
    normalEstimator.setParams(rayon / h);
    normalEstimator.init(h, surfels.begin(), surfels.end());

    for (int i = 0; i < surfels.size(); i++) {
        if (map(surfels[i]) == p) {
            RealPoint n = normalEstimator.eval(&(surfels[i]));
            cout << "Surfel : " << surfels[i] << " ; Normale : " << n << "\n";
            viewer_sphere.addBall(n, 0.0045);
            viewer_sphere.addBall(-n, 0.0045);
        }
    }
    cout << "\n";


    // estimation de la normale par plane probing
    cout << "Estimation de la normale par plane probing (mauve) :\n";
    viewer_sphere << CustomColors3D(Color(255, 128, 255), Color(255, 128, 255));

    using Surface = SH3::DigitalSurface;
    using Surfel = SH3::Surfel;
    DigitalSurfacePredicate<Surface> surf_pred = DigitalSurfacePredicate<Surface>(surface);

    Integer bound = params["maxAABB"].as<Integer>();
    double gridstep = params["gridstep"].as<double>();

    using ProbingAlgorithm = PlaneProbingParallelepipedEstimator<DigitalSurfacePredicate<Surface>, ProbingMode::R1>;
    using Estimator = PlaneProbingDigitalSurfaceLocalEstimator<Surface, ProbingAlgorithm>;

    Estimator::ProbingFactory probingFactory = [&bound](const Estimator::ProbingFrame& frame, const DigitalSurfacePredicate<Surface>& surfacePredicate) {
        return new ProbingAlgorithm(frame.p, { frame.b1, frame.b2, frame.normal }, surfacePredicate, bound);
    };

    std::unordered_map<Surfel, RealPoint> preEstimations;
    bool verbose = true;
    Estimator estimator(surface, probingFactory, preEstimations, verbose);
    estimator.init(gridstep, surfels.begin(), surfels.end());

    for (int i = 0; i < surfels.size(); i++) {
        if (map(surfels[i]) == p) {
            RealPoint n = estimator.eval(&surfels[i]).getNormalized() ;
            cout << "Surfel : " << surfels[i] << " ; Normale : " << n << "\n";
            viewer_sphere.addBall(n, 0.0035);
            viewer_sphere.addBall(-n, 0.0035);
        }
    }
    cout << "\n\n";


    // récupération de la normale réelle SUR LA SURFACE DE GOURSAT
    cout << "Normale réelle sur la surface de Goursat en " << p << " (blanc) :\n";
    viewer_sphere << CustomColors3D(Color::White, Color::White);
    params("polynomial", "goursat")("gridstep", h);
    auto implicit_shape = SH3::makeImplicitShape3D(params);
    auto normals = SHG3::getNormalVectors(implicit_shape, K, surfels, params);

    vector<RealPoint> nns;

    for (int i = 0; i < surfels.size(); i++) {
        if (map(surfels[i]) == p) {
            cout << "Surfel : " << surfels[i] << " ; Normale : " << normals[i] << "\n";
            viewer_sphere.addBall(normals[i], 0.005);
            viewer_sphere.addBall(-normals[i], 0.005);
            nns.push_back(normals[i]);
        }
    }
    viewer.addLine(p, p + 10 * BG(nns), 1);
    viewer << Point(p + 1);
            viewer << Viewer3D<>::updateDisplay;
    cout << "\n";


    // affichage des points de la grille de parcours considérés lors du calcul de BG Chord
    /*
    viewer_sphere << CustomColors3D(Color::Gray, Color::Gray);
    double pas = 0.01;
    RealPoint moyenne = RealPoint(0, 0, 0);
    double phi_min = std::numeric_limits<double>::max();

    double theta_min = 0;
    double theta_max = M_PI;
    double delta_min = - M_PI;
    double delta_max = M_PI;
    

    for (double theta = min(theta_min, theta_max); theta <= max(theta_min, theta_max); theta += pas) {

        double pas_delta = pas; 
        if (theta != 0) {
            pas_delta = pas / sin(theta);
        }
        for (double delta = min(delta_min, delta_max); delta <= max(delta_min, delta_max); delta += pas_delta) {

            RealPoint p = RealPoint(sin(theta) * cos(delta), sin(theta) * sin(delta), cos(theta)).getNormalized();
            double phi = 0;

            for (int i = 0; i < normales_Chord_decoups_BG.size(); i++) {
                double dot_prod = p.dot(normales_Chord_decoups_BG[i]);
                if (dot_prod < 1) {
                    phi += acos(dot_prod) * acos(dot_prod);
                }
            }
            if (phi < phi_min) {
                phi_min = phi;
                moyenne = p;
                viewer_sphere.addBall(p, 0.0025);
            }
        }
    }
    */

    // affichage des autres voxels de l'intersection
    // actuellement inutile tant qu'on affiche le dernier découpage
    /*
    for (auto s : surfels) {
        double dist = SP.distance(point_to_id[map(s)]);
        if (voxels_intersection(map(s))) {
            //cout << map(s) << "\n";
            //viewer << CustomColors3D(gradient(dist), gradient(dist));
            //viewer.setFillTransparency(207);
            //viewer << CustomColors3D(Color::Gray, Color::Gray);

            viewer << CustomColors3D(Color::Cyan, Color::Cyan);
            viewer.setFillTransparency(207);
        } else {
            viewer << CustomColors3D(Color::White, Color::White);
            //viewer.setFillTransparency(224);
        }
        viewer << map(s);
    }
    viewer << Viewer3D<>::updateDisplay;
    */

    viewer_sphere << Viewer3D<>::updateDisplay;
    application.exec();

    return 0;
}
