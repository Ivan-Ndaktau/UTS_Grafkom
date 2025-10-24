'use strict';
function createShiftryEllipsoid(x, y, z, rx, ry, rz, latBands, lonBands) {
    const vertices = [], normals = [], indices = [];
    for (let lat = 0; lat <= latBands; lat++) {
        const theta = lat * Math.PI / latBands;
        const sinTheta = Math.sin(theta);
        const cosTheta = Math.cos(theta);
        for (let lon = 0; lon <= lonBands; lon++) {
            const phi = lon * 2 * Math.PI / lonBands;
            const sinPhi = Math.sin(phi);
            const cosPhi = Math.cos(phi);
            const px = x + rx * cosPhi * sinTheta;
            const py = y + ry * cosTheta;
            const pz = z + rz * sinPhi * sinTheta;
            vertices.push(px, py, pz);
            const normal = VEC3.normalize([
                (px - x) / (rx*rx || 1e-6),
                (py - y) / (ry*ry || 1e-6),
                (pz - z) / (rz*rz || 1e-6)
            ]);
             if (isNaN(normal[0]) || VEC3.dot(normal,normal) < 0.1) {
                 const fallback = VEC3.normalize([px - x, py - y, pz - z]);
                 if(isNaN(fallback[0]) || VEC3.dot(fallback,fallback) < 0.1) {
                      normals.push(0, (py-y > 0 ? 1 : -1), 0);
                 } else {
                    normals.push(fallback[0], fallback[1], fallback[2]);
                 }
             } else {
                normals.push(normal[0], normal[1], normal[2]);
             }
        }
    }
    for (let lat = 0; lat < latBands; lat++) {
        for (let lon = 0; lon < lonBands; lon++) {
            const first = (lat * (lonBands + 1)) + lon;
            const second = first + lonBands + 1;
            indices.push(first, second, first + 1); // Winding CCW
            indices.push(second, second + 1, first + 1); // Winding CCW
        }
    }
    return { vertices, normals, indices };
}

// Membuat kerucut standar dengan pivot di tengah dasar (Y=0), cap tertutup, normal benar.
function createNuzleafCone(radius, height, segments = 12) {
    const vertices = [], normals = [], indices = [];
    const apex = [0, height, 0];       // Puncak di Y=height (relatif pivot)
    const baseCenter = [0, 0, 0];    // Dasar di Y=0 (relatif pivot)
    const angleStep = (2 * Math.PI) / segments;

    // 1. Apex vertex
    const apexIndex = 0;
    vertices.push(...apex);
    normals.push(0, 1, 0); // Normal puncak lurus ke atas

    // 2. Base vertices UNTUK SISI CONE (Normal miring)
    const baseVerticesSideStart = 1;
    for (let i = 0; i <= segments; i++) {
        const angle = i * angleStep;
        const cosA = Math.cos(angle);
        const sinA = Math.sin(angle);
        const baseX = cosA * radius;
        const baseZ = sinA * radius;
        vertices.push(baseX, 0, baseZ); // Titik dasar di Y=0

        // Normal Sisi Cone:
        // Vektor dari apex ke titik dasar
        const slopeVec = [baseX - apex[0], 0 - apex[1], baseZ - apex[2]]; // Y dasar = 0
        // Vektor tangen di dasar (tegak lurus radius)
        const tangentVec = [-sinA, 0, cosA]; // Arah CCW
        // Normal = cross(tangent, slope) lalu normalize & arahkan keluar
        let normal = VEC3.normalize(VEC3.cross(tangentVec, slopeVec));
        // Fallback jika normal tidak valid (misal radius=0 atau height=0)
        if(isNaN(normal[0]) || VEC3.dot(normal, normal) < 0.1) {
             normal = VEC3.normalize([cosA, Math.abs(radius/(height||1e-6)), sinA]);
        }
        normals.push(normal[0], normal[1], normal[2]);
    }

    // 3. Base vertices UNTUK CAP DASAR (Normal ke bawah)
    const baseCapCenterIndex = vertices.length / 3;
    vertices.push(...baseCenter);
    normals.push(0, -1, 0); // Normal tengah dasar ke bawah
    const baseCapRingStart = vertices.length / 3;
    for (let i = 0; i <= segments; i++) {
        const angle = i * angleStep;
        // Tambahkan vertex lagi untuk cap dengan normal berbeda
        vertices.push(Math.cos(angle) * radius, 0, Math.sin(angle) * radius);
        normals.push(0, -1, 0); // Normal ring dasar ke bawah
    }

    // 4. Indices
    for (let i = 0; i < segments; i++) {
        // Sisi Cone: Puncak, Titik Dasar Sisi (i+1), Titik Dasar Sisi (i) -> CCW dari luar
        indices.push(apexIndex, baseVerticesSideStart + i + 1, baseVerticesSideStart + i);
        // Cap Dasar: Tengah Dasar, Titik Cap Dasar (i), Titik Cap Dasar (i+1) -> CCW dari bawah
        indices.push(baseCapCenterIndex, baseCapRingStart + i, baseCapRingStart + i + 1);
    }

    return { vertices, normals, indices };
}

// Membuat silinder standar dengan pivot di tengah (Y=0), cap tertutup, normal benar.
function createNuzleafCylinder(radius, height, segments = 12) {
    const vertices = [], normals = [], indices = [];
    const halfHeight = height / 2;
    const angleStep = (2 * Math.PI) / segments;

    // Center vertices for caps
    const topCenterIndex = 0; vertices.push(0, halfHeight, 0); normals.push(0, 1, 0);
    const bottomCenterIndex = 1; vertices.push(0, -halfHeight, 0); normals.push(0, -1, 0);

    // Ring vertices (4 rings: topCap, bottomCap, sideTop, sideBottom)
    let currentVertIndex = 2; // Start after center vertices
    const topCapRingStart = currentVertIndex;
    for(let i=0; i<=segments; ++i){ vertices.push(radius * Math.cos(i*angleStep), halfHeight, radius * Math.sin(i*angleStep)); normals.push(0,1,0); }
    currentVertIndex += segments + 1;

    const bottomCapRingStart = currentVertIndex;
    for(let i=0; i<=segments; ++i){ vertices.push(radius * Math.cos(i*angleStep), -halfHeight, radius * Math.sin(i*angleStep)); normals.push(0,-1,0); }
    currentVertIndex += segments + 1;

    const sideRingStart = currentVertIndex;
    for(let i=0; i<=segments; ++i){
        const angle = i * angleStep;
        const cosA = Math.cos(angle); const sinA = Math.sin(angle);
        // Side Top vertex
        vertices.push(cosA * radius, halfHeight, sinA * radius); normals.push(cosA, 0, sinA);
        // Side Bottom vertex
        vertices.push(cosA * radius, -halfHeight, sinA * radius); normals.push(cosA, 0, sinA);
    }

    // Indices
    for (let i = 0; i < segments; i++) {
        // Top Cap
        indices.push(topCenterIndex, topCapRingStart + i, topCapRingStart + i + 1);
        // Bottom Cap (flipped winding)
        indices.push(bottomCenterIndex, bottomCapRingStart + i + 1, bottomCapRingStart + i);
        // Sides
        const sTop1 = sideRingStart + i*2;
        const sBot1 = sideRingStart + i*2 + 1;
        const sTop2 = sideRingStart + (i+1)*2;
        const sBot2 = sideRingStart + (i+1)*2 + 1;
        indices.push(sTop1, sBot1, sTop2); // CCW
        indices.push(sBot1, sBot2, sTop2); // CCW
    }
    return { vertices, normals, indices };
}

function createShiftryBox(x, y, z, width, height, depth) {
    const w=width/2, h=height/2, d=depth/2;
    const boxVertices=[ x-w,y-h,z+d, x+w,y-h,z+d, x+w,y+h,z+d, x-w,y+h,z+d, x+w,y-h,z-d, x-w,y-h,z-d, x-w,y+h,z-d, x+w,y+h,z-d, x-w,y+h,z+d, x+w,y+h,z+d, x+w,y+h,z-d, x-w,y+h,z-d, x-w,y-h,z-d, x+w,y-h,z-d, x+w,y-h,z+d, x-w,y-h,z+d, x+w,y-h,z+d, x+w,y-h,z-d, x+w,y+h,z-d, x+w,y+h,z+d, x-w,y-h,z-d, x-w,y-h,z+d, x-w,y+h,z+d, x-w,y+h,z-d ];
    const boxNormals=[ 0,0,1,0,0,1,0,0,1,0,0,1, 0,0,-1,0,0,-1,0,0,-1,0,0,-1, 0,1,0,0,1,0,0,1,0,0,1,0, 0,-1,0,0,-1,0,0,-1,0,0,-1,0, 1,0,0,1,0,0,1,0,0,1,0,0, -1,0,0,-1,0,0,-1,0,0,-1,0,0 ];
    const boxIndices = [ 0,1,2,0,2,3, 4,5,6,4,6,7, 8,9,10,8,10,11, 12,13,14,12,14,15, 16,17,18,16,18,19, 20,21,22,20,22,23 ];
    return { vertices: boxVertices, normals: boxNormals, indices: boxIndices };
}

// Tapered Cylinder dengan pivot di tengah Y=0, cap otomatis
function createNuzleafTaperedCylinder(topRadius, bottomRadius, height, segments) {
    const halfHeight = height / 2;
    const profile = [ { r: bottomRadius, y: -halfHeight }, { r: topRadius, y: halfHeight } ];
    const taperGeom = createGeometryFromProfile(profile, segments);
    const topCapCenterIndex = taperGeom.vertices.length / 3;
    taperGeom.vertices.push(0, halfHeight, 0); taperGeom.normals.push(0, 1, 0);
    const lastRingStartIndexTop = 0;
    for (let i = 0; i < segments; i++) {
        const tl = lastRingStartIndexTop + i * 2 + 1; const tr = lastRingStartIndexTop + (i + 1) * 2 + 1;
        taperGeom.indices.push(tl, topCapCenterIndex, tr);
    }
    const bottomCapCenterIndex = taperGeom.vertices.length / 3;
    taperGeom.vertices.push(0, -halfHeight, 0); taperGeom.normals.push(0, -1, 0);
    const firstRingStartIndexBottom = 0;
    for (let i = 0; i < segments; i++) {
        const bl = firstRingStartIndexBottom + i * 2; const br = firstRingStartIndexBottom + (i + 1) * 2;
        taperGeom.indices.push(br, bottomCapCenterIndex, bl);
    }
    return taperGeom;
}

function createNuzleafFoot() {
    const footLength = 0.4, footWidth = 0.25, footHeight = 0.15;
    const g1 = createShiftryBox(0, 0, 0, footWidth, footHeight, footLength);
    const g2 = createShiftryEllipsoid(0, 0, footLength/2, footWidth/2, footHeight/2, footHeight/2, 8, 8);
    const sideBlockSize = 0.08;
    const g3 = createShiftryBox(-footWidth/2 - sideBlockSize/2, 0, 0, sideBlockSize, footHeight, footLength/2);
    const g4 = createShiftryBox(footWidth/2 + sideBlockSize/2, 0, 0, sideBlockSize, footHeight, footLength/2);
    const vertices = [...g1.vertices, ...g2.vertices, ...g3.vertices, ...g4.vertices];
    const normals = [...g1.normals, ...g2.normals, ...g3.normals, ...g4.normals];
    const indices = [...g1.indices]; let offset = g1.vertices.length / 3;
    g2.indices.forEach(i => indices.push(i + offset)); offset += g2.vertices.length / 3;
    g3.indices.forEach(i => indices.push(i + offset)); offset += g3.vertices.length / 3;
    g4.indices.forEach(i => indices.push(i + offset));
    return { vertices, normals, indices };
}

function createNuzleafFlatShape(type, width, height, segments) {
    const vertices = [], normals = [], indices = []; const zOffset = 0;
    if (type === 'circle') { const r = width; vertices.push(0,0,zOffset); normals.push(0,0,1); for(let i=0;i<=segments;i++){ const a=(i/segments)*Math.PI*2; vertices.push(Math.cos(a)*r,Math.sin(a)*r,zOffset); normals.push(0,0,1); } for(let i=1;i<=segments;i++){ indices.push(0,i,i+1); } }
    else if (type === 'rightTriangle') { vertices.push(-width/2,-height/2,zOffset, width/2,-height/2,zOffset, -width/2,height/2,zOffset); normals.push(0,0,1,0,0,1,0,0,1); indices.push(0,1,2); }
    else if (type === 'rightTriangleMirrored') { vertices.push(-width/2,-height/2,zOffset, width/2,-height/2,zOffset, width/2,height/2,zOffset); normals.push(0,0,1,0,0,1,0,0,1); indices.push(0,1,2); }
    else if (type === 'triangle') { vertices.push(0,height/2,zOffset, -width/2,-height/2,zOffset, width/2,-height/2,zOffset); normals.push(0,0,1,0,0,1,0,0,1); indices.push(0,1,2); }
    else if (type === 'bezierLeaf') { const p0=[0,0,0], p1=[0,height*0.3,0.1], p2=[0,height*0.7,-0.1], p3=[0,height,0]; for(let i=0; i<=segments; i++) { const t=i/segments; const center=getBezierPoint(t,p0,p1,p2,p3); const tangent=getBezierTangent(t,p0,p1,p2,p3); const currentWidth=width*(1-Math.pow(Math.abs(2*t-1),2)); let binormal=Math.abs(tangent[1])<0.99?VEC3.normalize(VEC3.cross(tangent,[0,1,0])):[1,0,0]; let normal=VEC3.normalize(VEC3.cross(binormal,tangent)); if(isNaN(normal[0])) normal=[0,0,1]; const lV=VEC3.add(center,VEC3.scale(binormal,-currentWidth/2)); const rV=VEC3.add(center,VEC3.scale(binormal,currentWidth/2)); vertices.push(lV[0],lV[1],center[2]+zOffset); vertices.push(rV[0],rV[1],center[2]+zOffset); normals.push(...normal,...normal); } for(let i=0; i<segments; i++){ const i2=i*2; indices.push(i2,i2+2,i2+1); indices.push(i2+1,i2+2,i2+3); } }
    else if (type === 'strip') { const r=width, th=height; for(let i=0; i<=segments; i++){ const a=(i/segments)*Math.PI*2; const x=r*Math.cos(a), z=r*Math.sin(a); const xI=(r-th)*Math.cos(a), zI=(r-th)*Math.sin(a); vertices.push(x,0+zOffset,z, xI,0+zOffset,zI); normals.push(0,1,0, 0,1,0); } for(let i=0; i<segments; i++){ const i2=i*2; indices.push(i2,i2+2,i2+1); indices.push(i2+1,i2+2,i2+3); } }
    return { vertices, normals, indices };
}

// --- FUNGSI UTAMA PEMBUATAN GEOMETRI NUZLEAF ---
function createNuzleafGeometryData() {
    const SEGMENTS = 32;
    const headRadius = 0.7;
    const handRadius = 0.1;
    const handHeight = 0.25;

    return {
        head: createShiftryEllipsoid(0,0,0, headRadius*1.0, headRadius*0.9, headRadius*0.8, SEGMENTS, SEGMENTS),
        body: createNuzleafTaperedCylinder(0.35, 0.5, 0.8, SEGMENTS),
        nose: createNuzleafCone(0.1, 0.3, 20), // radius, height, segments (Pivot Y=0 di dasar)
        arm: createNuzleafCylinder(0.08, 1.2, 20), // radius, height, segments
        leg: createNuzleafCylinder(0.1, 1.8, 20),  // radius, height, segments
        leafStem: createNuzleafCylinder(0.04, 0.4, 10), // radius, height, segments
        hand: createShiftryEllipsoid(0,0,0, 0.1, 0.18, 0.15, 12, 12),
        diaperHalf: createShiftryEllipsoid(0,0,0, 0.5, 0.5, 0.5, 20, 20),
        mouth: createShiftryEllipsoid(0,0,0, 0.05, 0.05, 0.01, 10, 10),
        foot: createNuzleafFoot(),
        leftEye: createNuzleafFlatShape('rightTriangle', 0.3, 0.2, 0),
        rightEye: createNuzleafFlatShape('rightTriangleMirrored', 0.3, 0.2, 0),
        eyePupil: createNuzleafFlatShape('circle', 0.05, 0, 16),
        eyePupilWhite: createNuzleafFlatShape('circle', 0.015, 0, 12),
        diaperStrip: createNuzleafFlatShape('strip', 0.48, 0.02, 32),
        headStrip: createNuzleafFlatShape('strip', 0.62, 0.04, 32),
        leafBlade: createNuzleafFlatShape('bezierLeaf', 0.5, 1.2, SEGMENTS)
    };
}