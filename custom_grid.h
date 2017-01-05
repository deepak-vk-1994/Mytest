void populateFromCustom(std::vector<Point > &points, std::vector<Element > &elements) {
    ifstream file(gridfile);
    double coord1,coord2,coord3;

    Point pt;
    Element elem;
    int ne,np,nf;
    // //coarse
    if (gridfile == "./Grids/coarsecylinder.txt") {
        ne = 135;
        np = 216;
    }

    // //medium
    else if (gridfile == "./Grids/mediumcylinder.txt") {
        ne = 519;
        np = 808;
    }

    //fine
    else if (gridfile == "./Grids/finecylinder.txt") {
        ne = 2055;
        np = 3144;
    }

    //vfine
    else if (gridfile == "./Grids/vfinecylinder.txt") {
        ne = 8199;
        np = 12424;
    }

    else cout << "WHAT?";
    
    string line;
    int check = 1;
    double coord[3];
    int ptindex[4];
    getline(file,line);
    getline(file,line);
    getline(file,line);
    getline(file,line);
    getline(file,line);
    getline(file,line);
    getline(file,line);

    int index = 7;
    int rowindex = 0;
    while (index < ne) {
        getline(file,line);
        istringstream iss(line); 
        string word;
        rowindex = 0;
        while(iss >> word) {
            if (rowindex > 3) break;    
            check = 1;
            try {
                stod(word);   //string to double
            }
            catch(const std::invalid_argument& e) {
                check = 0;
            }
            if (check != 0) {
                ptindex[rowindex] = stod(word)-1;
                rowindex++;
            }
        }
        elem.vertex[0] = ptindex[1];
        elem.vertex[1] = ptindex[2];
        elem.vertex[2] = ptindex[3];
        elements.push_back(elem);
        index++;
    }
    getline(file,line);
    index++;
    while (index < np) {
        getline(file,line);
        istringstream iss(line); 
        string word;
        rowindex = 0;
        while(iss >> word) {
            if (rowindex > 2) break;    
            check = 1;
            try {
                stod(word);   //string to double
            }
            catch(const std::invalid_argument& e) {
                check = 0;
            }
            if (check != 0) {
                coord[rowindex] = stod(word);
                rowindex++;
            }
        }
        pt.x = coord[1]; 
        pt.y = coord[2];
        pt.z = 0.0;
        points.push_back(pt);
        index++;
    }

    file.close();
    N_elem = elements.size();
    N_pts = points.size();
    if (debug == 1) {
        ofstream ofile("./Debug/elements.txt");
        for (int i = 0; i < N_elem; i++) {
            // ofile << points[elements[i].vertex[0]].x<<" "<<points[elements[i].vertex[0]].y<<" "<<points[elements[i].vertex[0]].z<<endl;
            // ofile << points[elements[i].vertex[1]].x<<" "<<points[elements[i].vertex[1]].y<<" "<<points[elements[i].vertex[1]].z<<endl;
            // ofile << points[elements[i].vertex[2]].x<<" "<<points[elements[i].vertex[2]].y<<" "<<points[elements[i].vertex[2]].z<<endl;
            ofile << elements[i].vertex[0] << " " << elements[i].vertex[1] << " " << elements[i].vertex[2] << endl;
        }
        ofstream ofile2("./Debug/points.txt");
        for (int i = 0; i < N_pts; i++) {
            ofile2 << points[i].x<<" "<<points[i].y<<" "<<points[i].z<<endl;
        }
    }
    //ofile.close();
}
