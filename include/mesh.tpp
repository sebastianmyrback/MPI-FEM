







mesh1d::mesh1d(double a, double b, int n) {

    // n = number of vertices
    // n-1 = number of elements

    this->nv = n+1;
    this->nk = n;

    h = (b - a) / (nk);

    this->mesh_vertices.clear();
    this->elements.clear();
    this->mesh_vertices.resize(nv);
    this->elements.resize(nk);

    // Create vertices
    for (int i = 0; i < nv; i++) {
        const Rn x(a + i*h);  // x coordinate
        (this->mesh_vertices)[i].x            = x;
        (this->mesh_vertices)[i].glb_idx      = i;
        (this->mesh_vertices)[i].vertex_label = 0;
    }

    // Mark the boundary vertices
    (this->mesh_vertices)[0].vertex_label    = 1;
    (this->mesh_vertices)[nv-1].vertex_label = 2;

    border_dofs.push_back(0);
    border_dofs.push_back(nv-1);

    // Create elements
    for (int i = 0; i < nk; i++) {
        (this->elements)[i].elem_vertices = {std::make_shared<vert>(mesh_vertices[i]), std::make_shared<vert>(mesh_vertices[i + 1])};
        (this->elements)[i].measure = h;
    }



    return;
};