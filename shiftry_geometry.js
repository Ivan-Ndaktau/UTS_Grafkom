var VEC3_LIBS = {
    create: function (x = 0, y = 0, z = 0) { return [x, y, z]; },
    subtract: function (a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; },
    add: function (a, b) { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; },
    scale: function (v, s) { return [v[0] * s, v[1] * s, v[2] * s]; },
    normalize: function (v) {
        var len = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (len > 0.00001) { return [v[0] / len, v[1] / len, v[2] / len]; }
        else { return [0, 0, 0]; }
    },
    cross: function (a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; },
    dot: function (a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
};


var LIBS = {
    degToRad: function (angle) {
        return (angle * Math.PI / 180);
    },

    get_projection: function (angle, a, zMin, zMax) {
        var tan = Math.tan(LIBS.degToRad(0.5 * angle)),
            A = -(zMax + zMin) / (zMax - zMin),
            B = (-2 * zMax * zMin) / (zMax - zMin);

        return [
            0.5 / tan, 0, 0, 0,
            0, 0.5 * a / tan, 0, 0,
            0, 0, A, -1,
            0, 0, B, 0
        ];
    },

    get_I4: function () {
        return [1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1];
    },

    set_I4: function (m) {
        m[0] = 1, m[1] = 0, m[2] = 0, m[3] = 0,
            m[4] = 0, m[5] = 1, m[6] = 0, m[7] = 0,
            m[8] = 0, m[9] = 0, m[10] = 1, m[11] = 0,
            m[12] = 0, m[13] = 0, m[14] = 0, m[15] = 1;
    },

    // --- [BARU] Fungsi LookAt untuk Kamera FPS ---
    lookAt: function (eye, center, up) {
        var z = VEC3_LIBS.normalize(VEC3_LIBS.subtract(eye, center));
        var x = VEC3_LIBS.normalize(VEC3_LIBS.cross(up, z));
        var y = VEC3_LIBS.normalize(VEC3_LIBS.cross(z, x));

        return [
            x[0], y[0], z[0], 0,
            x[1], y[1], z[1], 0,
            x[2], y[2], z[2], 0,
            -VEC3_LIBS.dot(x, eye), -VEC3_LIBS.dot(y, eye), -VEC3_LIBS.dot(z, eye), 1
        ];
    },


    rotate: function (m, angle, axis) {
        var v = VEC3_LIBS.normalize(axis);
        var x = v[0], y = v[1], z = v[2];
        var c = Math.cos(angle);
        var s = Math.sin(angle);
        var C = 1 - c;

        var rot = [
            x * x * C + c, y * x * C + z * s, z * x * C - y * s, 0,
            x * y * C - z * s, y * y * C + c, z * y * C + x * s, 0,
            x * z * C + y * s, y * z * C - x * s, z * z * C + c, 0,
            0, 0, 0, 1
        ];

        var res = LIBS.multiply(m, rot);
        // Salin nilai kembali ke m
        for (var i = 0; i < 16; i++) {
            m[i] = res[i];
        }
    },

    rotateX: function (m, angle) {
        var c = Math.cos(angle);
        var s = Math.sin(angle);
        var mv1 = m[1], mv5 = m[5], mv9 = m[9];
        m[1] = m[1] * c - m[2] * s;
        m[5] = m[5] * c - m[6] * s;
        m[9] = m[9] * c - m[10] * s;

        m[2] = m[2] * c + mv1 * s;
        m[6] = m[6] * c + mv5 * s;
        m[10] = m[10] * c + mv9 * s;
    },

    rotateY: function (m, angle) {
        var c = Math.cos(angle);
        var s = Math.sin(angle);
        var mv0 = m[0], mv4 = m[4], mv8 = m[8];
        m[0] = c * m[0] + s * m[2];
        m[4] = c * m[4] + s * m[6];
        m[8] = c * m[8] + s * m[10];

        m[2] = c * m[2] - s * mv0;
        m[6] = c * m[6] - s * mv4;
        m[10] = c * m[10] - s * mv8;
    },

    rotateZ: function (m, angle) {
        var c = Math.cos(angle);
        var s = Math.sin(angle);
        var mv0 = m[0], mv4 = m[4], mv8 = m[8];
        m[0] = c * m[0] - s * m[1];
        m[4] = c * m[4] - s * m[5];
        m[8] = c * m[8] - s * m[9];

        m[1] = c * m[1] + s * mv0;
        m[5] = c * m[5] + s * mv4;
        m[9] = c * m[9] + s * mv8;
    },

    translateZ: function (m, t) {
        m[14] += t;
    },

    translateX: function (m, t) {
        m[12] += t;
    },

    translateY: function (m, t) {
        m[13] += t;
    },

    set_position: function (m, x, y, z) {
        m[12] = x, m[13] = y, m[14] = z;
    },

    scaleX: function (m, t) {
        m[0] *= t;
    },

    scaleY: function (m, t) {
        m[5] *= t;
    },

    scaleZ: function (m, t) {
        m[10] *= t;
    },

    multiply: function (m1, m2) {
        var rm = this.get_I4();
        var N = 4;
        for (var i = 0; i < N; i++) {
            for (var j = 0; j < N; j++) {
                rm[i * N + j] = 0;
                for (var k = 0; k < N; k++)
                    rm[i * N + j] += m1[i * N + k] * m2[k * N + j];
            }
        }
        return rm;
    },

    scale: function (m, s) {
        if (Array.isArray(s)) {
            m[0] *= s[0];
            m[5] *= s[1];
            m[10] *= s[2];
        } else {
            m[0] *= s;
            m[5] *= s;
            m[10] *= s;
        }
    },

    copy: function (src) {
        var dst = [];
        for (var i = 0; i < 16; i++) {
            dst[i] = src[i];
        }
        return dst;
    }
};



var VEC3 = VEC3_LIBS;

function matrixInverse(m) {
    var inv = new Array(16);
    inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
    inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
    inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
    inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
    inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
    inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
    inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
    inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
    inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
    inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
    inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
    inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
    inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
    inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
    inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
    inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

    var det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
    if (Math.abs(det) < 1e-15) return LIBS.get_I4();
    det = 1.0 / det;
    for (var i = 0; i < 16; i++) inv[i] = inv[i] * det;
    return inv;
}
function matrixTranspose(m) { return [m[0], m[4], m[8], m[12], m[1], m[5], m[9], m[13], m[2], m[6], m[10], m[14], m[3], m[7], m[11], m[15]]; }
function get_Nmatrix(m) { return matrixTranspose(matrixInverse(m)); }

VEC3.transformMat4Normal = function (v, m) {
    var x = v[0], y = v[1], z = v[2];
    return [
        m[0] * x + m[4] * y + m[8] * z,
        m[1] * x + m[5] * y + m[9] * z,
        m[2] * x + m[6] * y + m[10] * z
    ];
};
VEC3.transformMat4Position = function (v, m) {
    var x = v[0], y = v[1], z = v[2];
    var w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
    return [
        (m[0] * x + m[4] * y + m[8] * z + m[12]) / w,
        (m[1] * x + m[5] * y + m[9] * z + m[13]) / w,
        (m[2] * x + m[6] * y + m[10] * z + m[14]) / w
    ];
};


function getBezierPoint(t, p0, p1, p2, p3) {
    var cX = 3 * (p1[0] - p0[0]), bX = 3 * (p2[0] - p1[0]) - cX, aX = p3[0] - p0[0] - cX - bX;
    var cY = 3 * (p1[1] - p0[1]), bY = 3 * (p2[1] - p1[1]) - cY, aY = p3[1] - p0[1] - cY - bY;
    var cZ = 3 * (p1[2] - p0[2]), bZ = 3 * (p2[2] - p1[2]) - cZ, aZ = p3[2] - p0[2] - cZ - bZ;
    var t2 = t * t, t3 = t2 * t;
    return [(aX * t3) + (bX * t2) + (cX * t) + p0[0],
    (aY * t3) + (bY * t2) + (cY * t) + p0[1],
    (aZ * t3) + (bZ * t2) + (cZ * t) + p0[2]];
}
function getBezierTangent(t, p0, p1, p2, p3) {
    var cX = 3 * (p1[0] - p0[0]), bX = 3 * (p2[0] - p1[0]) - cX, aX = p3[0] - p0[0] - cX - bX;
    var cY = 3 * (p1[1] - p0[1]), bY = 3 * (p2[1] - p1[1]) - cY, aY = p3[1] - p0[1] - cY - bY;
    var cZ = 3 * (p1[2] - p0[2]), bZ = 3 * (p2[2] - p1[2]) - cZ, aZ = p3[2] - p0[2] - cZ - bZ;
    var t2 = t * t;
    return VEC3.normalize([(3 * aX * t2) + (2 * bX * t) + cX,
    (3 * aY * t2) + (2 * bY * t) + cY,
    (3 * aZ * t2) + (2 * bZ * t) + cZ]);
}
function createCurvedShape(controlPoints, segments, radius, color, taper = 0.9, radialSegments = 6, flatten = 1.0) {
    var vertices = [], indices = [], colors = [], normals = [];
    var lastBinormal = [1, 0, 0];
    for (var i = 0; i <= segments; i++) {
        var t = i / segments;
        var center = getBezierPoint(t, controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3]);
        var tangent = getBezierTangent(t, controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3]);
        var approxUp = [0, 1, 0];
        if (Math.abs(VEC3.dot(tangent, approxUp)) > 0.99) {
            approxUp = lastBinormal;
        }
        var binormal = VEC3.normalize(VEC3.cross(tangent, approxUp));
        var normal = VEC3.normalize(VEC3.cross(binormal, tangent));
        lastBinormal = binormal;
        var currentRadius = radius * (1.0 - (t * taper));
        for (var j = 0; j <= radialSegments; j++) {
            var angle = (j / radialSegments) * 2.0 * Math.PI;
            var cosA = Math.cos(angle);
            var sinA = Math.sin(angle);
            var offsetVec = VEC3.add(
                VEC3.scale(normal, cosA * currentRadius * flatten),
                VEC3.scale(binormal, sinA * currentRadius)
            );
            var vX = center[0] + offsetVec[0];
            var vY = center[1] + offsetVec[1];
            var vZ = center[2] + offsetVec[2];
            vertices.push(vX, vY, vZ);
            colors.push(color[0], color[1], color[2]);
            var vertexNormal = VEC3.normalize(offsetVec);
            if (vertexNormal[0] === 0 && vertexNormal[1] === 0 && vertexNormal[2] === 0) {
                vertexNormal = normal;
            }
            normals.push(vertexNormal[0], vertexNormal[1], vertexNormal[2]);
        }
    }
    for (var i = 0; i < segments; i++) {
        for (var j = 0; j < radialSegments; j++) {
            var first = i * (radialSegments + 1) + j;
            var second = first + radialSegments + 1;
            indices.push(first, second, first + 1);
            indices.push(second, second + 1, first + 1);
        }
    }
    var lastCenter = getBezierPoint(1.0, controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3]);
    var tipIndexStart = vertices.length / 3;
    vertices.push(lastCenter[0], lastCenter[1], lastCenter[2]);
    colors.push(color[0] * 0.9, color[1] * 0.9, color[2] * 0.9);
    normals.push(0, -1, 0);
    var baseIndexStart = (segments - 1) * (radialSegments + 1);
    for (var j = 0; j < radialSegments; j++) {
        indices.push(tipIndexStart, baseIndexStart + j, baseIndexStart + j + 1);
    }
    return { vertices: vertices, colors: colors, normals: normals, indices: indices };
}
function createEllipsoid(radiusX, radiusY, radiusZ, latitudeBands, longitudeBands, color) {
    var vertices = [], colors = [], indices = [], normals = [];
    for (var lat = 0; lat <= latitudeBands; lat++) {
        var theta = lat * Math.PI / latitudeBands;
        var sinTheta = Math.sin(theta); var cosTheta = Math.cos(theta);
        for (var lon = 0; lon <= longitudeBands; lon++) {
            var phi = lon * 2 * Math.PI / longitudeBands;
            var sinPhi = Math.sin(phi); var cosPhi = Math.cos(phi);
            var x = radiusX * cosPhi * sinTheta;
            var y = radiusY * cosTheta;
            var z = radiusZ * sinPhi * sinTheta;
            var norm = VEC3.normalize([x / (radiusX * radiusX), y / (radiusY * radiusY), z / (radiusZ * radiusZ)]);
            if (isNaN(norm[0]) || isNaN(norm[1]) || isNaN(norm[2])) norm = [0, y > 0 ? 1 : -1, 0];
            vertices.push(x, y, z);
            colors.push(color[0], color[1], color[2]);
            normals.push(norm[0], norm[1], norm[2]);
        }
    }
    for (var lat = 0; lat < latitudeBands; lat++) {
        for (var lon = 0; lon < longitudeBands; lon++) {
            var first = (lat * (longitudeBands + 1)) + lon;
            var second = first + longitudeBands + 1;
            indices.push(first, second, first + 1);
            indices.push(second, second + 1, first + 1);
        }
    }
    return { vertices: vertices, colors: colors, normals: normals, indices: indices };
}
function createSphere(radius, latitudeBands, longitudeBands, color) {
    return createEllipsoid(radius, radius, radius, latitudeBands, longitudeBands, color);
}
function createEllipticParaboloid(radiusA, radiusB, height, segmentsU, segmentsV, color, inverse = false) {
    var vertices = [], colors = [], indices = [], normals = [];
    var direction = inverse ? -1 : 1;
    for (var i = 0; i <= segmentsU; i++) {
        var u = i / segmentsU;
        var u_eps = 0.0001;
        var u_next = Math.min(u + u_eps, 1.0);
        for (var j = 0; j <= segmentsV; j++) {
            var v = j / segmentsV * 2 * Math.PI;
            var v_eps = 0.0001;
            var v_next = v + v_eps;
            var x = radiusA * u * Math.cos(v);
            var z = radiusB * u * Math.sin(v);
            var y = height * u * u * direction;
            vertices.push(x, y, z);
            colors.push(color[0], color[1], color[2]);
            var pt = [x, y, z];
            var pt_u = [radiusA * u_next * Math.cos(v), height * u_next * u_next * direction, radiusB * u_next * Math.sin(v)];
            var pt_v = [radiusA * u * Math.cos(v_next), y, radiusB * u * Math.sin(v_next)];
            var tangentU = VEC3.subtract(pt_u, pt);
            var tangentV = VEC3.subtract(pt_v, pt);
            var norm = VEC3.normalize(VEC3.cross(tangentU, tangentV));
            if (inverse) {
                norm = VEC3.scale(norm, -1.0);
            }
            var lenSq = norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2];
            if (lenSq < 0.00001) {
                norm = [0, direction, 0];
            }
            normals.push(norm[0], norm[1], norm[2]);
        }
    }
    for (var i = 0; i < segmentsU; i++) {
        for (var j = 0; j < segmentsV; j++) {
            var first = i * (segmentsV + 1) + j;
            var second = first + segmentsV + 1;
            indices.push(first, first + 1, second);
            indices.push(second, first + 1, second + 1);
        }
    }
    return { vertices: vertices, colors: colors, normals: normals, indices: indices };
}
function createDetailedLeaf(leafLength, leafWidth, leafThickness, color, segments = 12, widthSegments = 6) {
    var vertices = [], colors = [], normals = [], indices = [];
    function getWidth(t) { return leafWidth * Math.sin(t * Math.PI); }
    function getThickness(t) { return leafThickness * Math.sin(t * Math.PI); }
    function getBend(t) { return -leafWidth * 0.3 * Math.pow(Math.sin(t * Math.PI), 2); }
    for (var i = 0; i <= segments; i++) {
        var t = i / segments;
        var currentY = t * leafLength;
        var currentWidth = getWidth(t);
        var currentThickness = getThickness(t);
        var currentBend = getBend(t);
        for (var j = 0; j <= widthSegments; j++) {
            var s = j / widthSegments;
            var currentX = (s - 0.5) * currentWidth;
            var t_start = 0.00;
            var t_end = 0.5;
            var maxVeinWidth = 0.25;
            var veinColor = COLORS.BODY_BROWN;
            var currentVeinWidth = 0;
            if (t >= t_start && t <= t_end) {
                var taperFactor = 1.0 - ((t - t_start) / (t_end - t_start));
                currentVeinWidth = maxVeinWidth * taperFactor;
            }
            var isVeinArea = Math.abs(s - 0.5) < currentVeinWidth;
            var currentTopColor = color;
            var currentBottomColor = [color[0] * 0.8, color[1] * 0.8, color[2] * 0.8];
            var topZ = currentBend + currentThickness * 0.5 * Math.cos((s - 0.5) * Math.PI);
            var bottomZ = currentBend - currentThickness * 0.5 * Math.cos((s - 0.5) * Math.PI);
            if (isVeinArea) {
                currentTopColor = veinColor;
                currentBottomColor = veinColor;
            }
            vertices.push(currentX, currentY, topZ);
            colors.push(currentTopColor[0], currentTopColor[1], currentTopColor[2]);
            var normalTop = VEC3.normalize([(s - 0.5) * 0.2, -0.1, 1.0]);
            normals.push(normalTop[0], normalTop[1], normalTop[2]);
            vertices.push(currentX, currentY, bottomZ);
            colors.push(currentBottomColor[0], currentBottomColor[1], currentBottomColor[2]);
            var normalBottom = VEC3.normalize([(s - 0.5) * 0.2, 0.1, -1.0]);
            normals.push(normalBottom[0], normalBottom[1], normalBottom[2]);
        }
    }
    var pointsPerSegment = (widthSegments + 1) * 2;
    for (var i = 0; i < segments; i++) {
        for (var j = 0; j < widthSegments; j++) {
            var topLeft = i * pointsPerSegment + j * 2;
            var bottomLeft = topLeft + 1;
            var topRight = topLeft + 2;
            var bottomRight = topRight + 1;
            var nextTopLeft = (i + 1) * pointsPerSegment + j * 2;
            var nextBottomLeft = nextTopLeft + 1;
            var nextTopRight = nextTopLeft + 2;
            var nextBottomRight = nextTopRight + 1;
            indices.push(topLeft, nextTopLeft, topRight);
            indices.push(nextTopLeft, nextTopRight, topRight);
            indices.push(bottomLeft, bottomRight, nextBottomLeft);
            indices.push(nextBottomLeft, bottomRight, nextBottomRight);
            if (j == 0) {
                indices.push(bottomLeft, topLeft, nextBottomLeft);
                indices.push(topLeft, nextTopLeft, nextBottomLeft);
            }
            if (j == widthSegments - 1) {
                indices.push(bottomRight, nextBottomRight, topRight);
                indices.push(topRight, nextBottomRight, nextTopRight);
            }
        }
    }
    return { vertices: vertices, colors: colors, normals: normals, indices: indices };
}
function createCylinder(topRadius, bottomRadius, height, segments, topColor, bottomColor) {
    var vertices = [], colors = [], indices = [], normals = [];
    var topCenterIndex = 0; var bottomCenterIndex = 1;
    vertices.push(0, height / 2, 0); colors.push(topColor[0], topColor[1], topColor[2]); normals.push(0, 1, 0);
    vertices.push(0, -height / 2, 0); colors.push(bottomColor[0], bottomColor[1], bottomColor[2]); normals.push(0, -1, 0);
    var capVertStart = vertices.length / 3;
    for (var i = 0; i <= segments; i++) {
        var angle = (i / segments) * Math.PI * 2;
        var x = Math.cos(angle); var z = Math.sin(angle);
        vertices.push(x * topRadius, height / 2, z * topRadius); colors.push(topColor[0], topColor[1], topColor[2]); normals.push(0, 1, 0);
        vertices.push(x * bottomRadius, -height / 2, z * bottomRadius); colors.push(bottomColor[0], bottomColor[1], bottomColor[2]); normals.push(0, -1, 0);
    }
    var sideVertStart = vertices.length / 3;
    for (var i = 0; i <= segments; i++) {
        var angle = (i / segments) * Math.PI * 2;
        var x = Math.cos(angle); var z = Math.sin(angle);
        var deltaY = height; var deltaR = bottomRadius - topRadius;
        var norm = VEC3.normalize([x * deltaY, deltaR, z * deltaY]);
        vertices.push(x * topRadius, height / 2, z * topRadius); colors.push(topColor[0], topColor[1], topColor[2]); normals.push(norm[0], norm[1], norm[2]);
        vertices.push(x * bottomRadius, -height / 2, z * bottomRadius); colors.push(bottomColor[0], bottomColor[1], bottomColor[2]); normals.push(norm[0], norm[1], norm[2]);
    }
    for (var i = 0; i < segments; i++) {
        var top1 = capVertStart + i * 2; var top2 = capVertStart + (i + 1) * 2;
        var bot1 = capVertStart + i * 2 + 1; var bot2 = capVertStart + (i + 1) * 2 + 1;
        indices.push(topCenterIndex, top2, top1);
        indices.push(bottomCenterIndex, bot1, bot2);
    }
    for (var i = 0; i < segments; i++) {
        var sideTop1 = sideVertStart + i * 2; var sideBot1 = sideVertStart + i * 2 + 1;
        var sideTop2 = sideVertStart + (i + 1) * 2; var sideBot2 = sideVertStart + (i + 1) * 2 + 1;
        indices.push(sideTop1, sideBot1, sideTop2);
        indices.push(sideBot1, sideBot2, sideTop2);
    }
    return { vertices: vertices, colors: colors, normals: normals, indices: indices };
}
function createCone(radius, height, segments, color) {
    return createCylinder(0.001, radius, height, segments, color, color);
}
function createHyperboloid(waistRadius, topRadius, height, segmentsY, segmentsX, color) {
    var vertices = [], colors = [], indices = [], normals = [];
    var h2 = height / 2;
    var a = waistRadius;
    var c_squared;
    if (topRadius > waistRadius) {
        c_squared = (h2 * h2) / (Math.pow(topRadius / waistRadius, 2) - 1);
    } else {
        c_squared = 10000;
        console.warn("Hyperboloid topRadius <= waistRadius, approximation used.");
    }
    for (var i = 0; i <= segmentsY; i++) {
        var t = i / segmentsY;
        var y = height * t - h2;
        var radius = Math.sqrt(a * a * (1 + (y * y / c_squared)));
        for (var j = 0; j <= segmentsX; j++) {
            var s = j / segmentsX;
            var angle = s * 2 * Math.PI;
            var x = radius * Math.cos(angle);
            var z = radius * Math.sin(angle);
            vertices.push(x, y, z);
            colors.push(color[0], color[1], color[2]);
            var norm = VEC3.normalize([
                x / (a * a),
                -y / c_squared,
                z / (a * a)
            ]);
            if (isNaN(norm[0])) norm = VEC3.normalize([x, 0, z]);
            normals.push(norm[0], norm[1], norm[2]);
        }
    }
    for (var i = 0; i < segmentsY; i++) {
        for (var j = 0; j < segmentsX; j++) {
            var first = i * (segmentsX + 1) + j;
            var second = first + segmentsX + 1;
            indices.push(first, second, first + 1);
            indices.push(second, second + 1, first + 1);
        }
    }
    return { vertices: vertices, colors: colors, normals: normals, indices: indices };
}
function createBox(width, height, depth, color) {
    var w = width / 2, h = height / 2, d = depth / 2;
    var vertices = [-w, -h, d, w, -h, d, w, h, d, -w, h, d, -w, -h, -d, -w, h, -d, w, h, -d, w, -h, -d, -w, h, d, w, h, d, w, h, -d, -w, h, -d, -w, -h, d, -w, -h, -d, w, -h, -d, w, -h, d, w, -h, d, w, -h, -d, w, h, -d, w, h, d, -w, -h, d, -w, h, d, -w, h, -d, -w, -h, -d];
    var faceColors = [[color[0] * 1.0, color[1] * 1.0, color[2] * 1.0], [color[0] * 0.9, color[1] * 0.9, color[2] * 0.9], [color[0] * 0.8, color[1] * 0.8, color[2] * 0.8], [color[0] * 0.7, color[1] * 0.7, color[2] * 0.7], [color[0] * 0.6, color[1] * 0.6, color[2] * 0.6], [color[0] * 0.5, color[1] * 0.5, color[2] * 0.5]];
    var colors = []; for (var j = 0; j < faceColors.length; j++) { for (var i = 0; i < 4; i++) { colors = colors.concat(faceColors[j]); } }
    var normals = [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0];
    var indices = [0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 8, 9, 10, 8, 10, 11, 12, 13, 14, 12, 14, 15, 16, 17, 18, 16, 18, 19, 20, 21, 22, 20, 22, 23];
    return { vertices: vertices, colors: colors, normals: normals, indices: indices };
}

function createTorus(mainRadius, tubeRadius, mainSegments, tubeSegments, color) {
    var vertices = [], colors = [], normals = [], indices = [];

    for (var i = 0; i <= mainSegments; i++) {
        var u = i / mainSegments * 2 * Math.PI;
        var cosU = Math.cos(u);
        var sinU = Math.sin(u);

        var center = [mainRadius * cosU, 0, mainRadius * sinU];

        for (var j = 0; j <= tubeSegments; j++) {
            var v = j / tubeSegments * 2 * Math.PI;
            var cosV = Math.cos(v);
            var sinV = Math.sin(v);

            var x = (mainRadius + tubeRadius * cosV) * cosU;
            var y = tubeRadius * sinV;
            var z = (mainRadius + tubeRadius * cosV) * sinU;
            vertices.push(x, y, z);

            colors.push(color[0], color[1], color[2]);

            var tubeNormal = [
                tubeRadius * cosV * cosU,
                tubeRadius * sinV,
                tubeRadius * cosV * sinU
            ];
            var norm = VEC3.normalize(tubeNormal);
            if (isNaN(norm[0]) || VEC3.dot(norm, norm) < 0.001) norm = VEC3.normalize([cosU, sinV, sinU]);

            normals.push(norm[0], norm[1], norm[2]);
        }
    }

    for (var i = 0; i < mainSegments; i++) {
        for (var j = 0; j < tubeSegments; j++) {
            var first = (i * (tubeSegments + 1)) + j;
            var second = first + tubeSegments + 1;

            indices.push(first, second, first + 1);
            indices.push(second, second + 1, first + 1);
        }
    }

    return { vertices: vertices, colors: colors, normals: normals, indices: indices };
}

function combineGeometry(dataArray, transformArray) {
    var allVertices = [], allColors = [], allIndices = [], allNormals = [];
    var currentVertexOffset = 0;
    for (var i = 0; i < dataArray.length; i++) {
        var data = dataArray[i];
        if (!data || !data.vertices || !data.indices || !data.normals) {
            console.error("Data geometri tidak valid atau kehilangan normal/vertices/indices di indeks:", i, data); continue;
        }
        var transform = transformArray[i]; var partMatrix = LIBS.get_I4();
        LIBS.translateX(partMatrix, transform.x || 0); LIBS.translateY(partMatrix, transform.y || 0); LIBS.translateZ(partMatrix, transform.z || 0);
        if (transform.rotY) LIBS.rotateY(partMatrix, LIBS.degToRad(transform.rotY));
        if (transform.rotX) LIBS.rotateX(partMatrix, LIBS.degToRad(transform.rotX));
        if (transform.rotZ) LIBS.rotateZ(partMatrix, LIBS.degToRad(transform.rotZ));
        if (transform.scale) LIBS.scale(partMatrix, Array.isArray(transform.scale) ? transform.scale : [transform.scale, transform.scale, transform.scale]);

        var normalMatrix = get_Nmatrix(partMatrix);

        for (var j = 0; j < data.vertices.length; j += 3) {
            var v = [data.vertices[j], data.vertices[j + 1], data.vertices[j + 2]];
            var n = [data.normals[j], data.normals[j + 1], data.normals[j + 2]];

            var transformedV = VEC3.transformMat4Position(v, partMatrix);
            var transformedN = VEC3.normalize(VEC3.transformMat4Normal(n, normalMatrix));

            allVertices.push(transformedV[0], transformedV[1], transformedV[2]);
            allNormals.push(transformedN[0], transformedN[1], transformedN[2]);
        }

        var numVertices = data.vertices.length / 3;
        if (data.colors && data.colors.length === numVertices * 3) {
            allColors = allColors.concat(data.colors);
        } else {
            console.warn("Missing atau warna mismatch untuk geometry index:", i, ". Using white.");
            for (var k = 0; k < numVertices; ++k) allColors.push(1.0, 1.0, 1.0);
        }

        for (var j = 0; j < data.indices.length; j++) allIndices.push(data.indices[j] + currentVertexOffset);
        currentVertexOffset += numVertices;
    }
    return { vertices: allVertices, colors: allColors, normals: allNormals, indices: allIndices };
}

var COLORS = {
    BODY_BROWN: [0.55, 0.40, 0.32],
    HAIR_WHITE: [0.95, 0.95, 0.95],
    HAIR_SHADOW: [0.8, 0.8, 0.8],
    LEAF_GREEN: [0.2, 0.6, 0.2],
    LEAF_DARK: [0.1, 0.4, 0.1],
    EYE_YELLOW: [1.0, 0.9, 0.1],
    PUPIL_BLACK: [0.1, 0.1, 0.1],
    DARK_PULSE: [0.3, 0.0, 0.5],
    SWORD_BLUE: [0.3, 0.7, 1.0],
    SWORD_HILT: [0.8, 0.8, 0.8]
};

var torsoData = createEllipsoid(0.5, 0.6, 0.4, 30, 30, COLORS.BODY_BROWN);
var noseData = createCone(0.1, 0.6, 12, COLORS.BODY_BROWN);
var mouthData = createBox(0.3, 0.06, 0.02, COLORS.PUPIL_BLACK);
var toothData = createBox(0.3, 0.025, 0.02, COLORS.HAIR_WHITE);
var outerEarData = createEllipticParaboloid(0.15, 0.05, 0.6, 12, 12, COLORS.BODY_BROWN, true);
var innerEarData = createEllipticParaboloid(0.1, 0.05, 0.45, 10, 10, COLORS.PUPIL_BLACK, true);
var hairBaseData = createEllipsoid(0.5, 0.5, 0.4, 20, 20, COLORS.HAIR_WHITE);

function createCustomEye(centerX, centerY, mirrorX = false) {
    var eyeZBase = 0.38;
    var eyeLayerOffset = 0.0001;
    var pupilRadius = 0.02;
    var topY = 0.07;
    var bottomY = -0.07;
    var innerX_offset = 0.06;
    var outerX_offset = 0.15;
    var innerScale = 0.8;
    var outer_InnerX, outer_OuterX;
    var inner_InnerX, inner_OuterX;
    if (mirrorX) {
        outer_InnerX = centerX - innerX_offset;
        outer_OuterX = centerX + outerX_offset;
    } else {
        outer_InnerX = centerX + innerX_offset;
        outer_OuterX = centerX - outerX_offset;
    }
    var innerInnerX_offset = innerX_offset * innerScale;
    var innerOuterX_offset = outerX_offset * innerScale;
    if (mirrorX) {
        inner_InnerX = centerX - innerInnerX_offset;
        inner_OuterX = centerX + innerOuterX_offset;
    } else {
        inner_InnerX = centerX + innerInnerX_offset;
        inner_OuterX = centerX - innerOuterX_offset;
    }
    var outerVerts = [
        outer_InnerX, centerY + bottomY, eyeZBase,
        outer_OuterX, centerY + bottomY, eyeZBase,
        outer_OuterX, centerY + topY, eyeZBase
    ];
    var outerColors = [...COLORS.PUPIL_BLACK, ...COLORS.PUPIL_BLACK, ...COLORS.PUPIL_BLACK];
    var outerNormals = [0, 0, 1, 0, 0, 1, 0, 0, 1];
    var outerIndices = [0, 1, 2];
    var base_yellow_Y_offset = -0.015;
    var base_yellow_X_offset = 0.01;
    var base_pupil_Y_offset = 0.02;
    var base_pupil_X_offset = 0.09;
    var final_yellow_Y_offset = base_yellow_Y_offset;
    var final_pupil_Y_offset = base_pupil_Y_offset;
    var final_yellow_X_offset = mirrorX ? -base_yellow_X_offset : base_yellow_X_offset;
    var final_pupil_X_offset = mirrorX ? -base_pupil_X_offset : base_pupil_X_offset;
    var innerTopY_val = topY * innerScale;
    var innerBottomY_val = bottomY * innerScale;
    var innerVerts = [
        inner_InnerX + final_yellow_X_offset, centerY + innerBottomY_val + final_yellow_Y_offset, eyeZBase + eyeLayerOffset,
        inner_OuterX + final_yellow_X_offset, centerY + innerBottomY_val + final_yellow_Y_offset, eyeZBase + eyeLayerOffset,
        inner_OuterX + final_yellow_X_offset, centerY + innerTopY_val + final_yellow_Y_offset, eyeZBase + eyeLayerOffset
    ];
    var innerColors = [...COLORS.EYE_YELLOW, ...COLORS.EYE_YELLOW, ...COLORS.EYE_YELLOW];
    var innerNormals = [0, 0, 1, 0, 0, 1, 0, 0, 1];
    var innerIndices = [3, 4, 5];
    var flatEyeData = {
        vertices: [...outerVerts, ...innerVerts],
        colors: [...outerColors, ...innerColors],
        normals: [...outerNormals, ...innerNormals],
        indices: [...outerIndices, ...innerIndices]
    };
    var pupilHeight = 0.01;
    var pupilData = createCylinder(pupilRadius, pupilRadius, pupilHeight, 12, COLORS.PUPIL_BLACK, COLORS.PUPIL_BLACK);
    var pupilPosX = inner_OuterX + final_yellow_X_offset + final_pupil_X_offset;
    var pupilPosY = (centerY + innerBottomY_val * 0.85) + final_yellow_Y_offset + final_pupil_Y_offset;
    var pupilPosZ = eyeZBase + 2 * eyeLayerOffset;
    var finalEyeData = combineGeometry(
        [flatEyeData, pupilData],
        [
            { x: 0, y: 0, z: 0 },
            { x: pupilPosX, y: pupilPosY, z: pupilPosZ, rotX: 90 }
        ]
    );
    return finalEyeData;
}
var leftEyeData = createCustomEye(0.05, 0.25);
var rightEyeData = createCustomEye(-0.05, 0.25, true);

var hairShortCP = [[0, 0, 0], [0.0, -0.25, -0.3], [0.0, -0.6, -0.3], [-0.0, -1.2, -0.3]];
var hairShortData = createCurvedShape(hairShortCP, 10, 0.18, COLORS.HAIR_WHITE, 1.0, 6, 0.3);
var hairLongCP = [[0, 0, 0], [0.0, -0.5, -0.3], [0.0, -1.0, -0.3], [0.0, -1.5, -0.3]];
var hairLongData = createCurvedShape(hairLongCP, 12, 0.18, COLORS.HAIR_WHITE, 1.0, 6, 0.3);
var shoulderData = createSphere(0.12, 12, 12, COLORS.BODY_BROWN);
var armData = createCone(0.05, 0.3, 8, COLORS.BODY_BROWN);
var leafData = createDetailedLeaf(0.7, 0.30, 0.2, COLORS.LEAF_GREEN, 12, 100);
var leafDataMid = createDetailedLeaf(0.85, 0.30, 0.2, COLORS.LEAF_GREEN, 12, 100);
var thighData = createEllipsoid(0.2, 0.3, 0.25, 12, 12, COLORS.BODY_BROWN);
var legData = createHyperboloid(0.05, 0.1, 0.15, 10, 12, COLORS.BODY_BROWN);
var footData = createBox(0.2, 0.1, 0.45, COLORS.BODY_BROWN);
var heelData = createBox(0.2, 0.15, 0.1, COLORS.BODY_BROWN);

var hairClumpData = [];
var hairClumpTransforms = [];
hairClumpData.push(
    hairShortData, hairShortData, hairShortData, hairShortData, hairShortData,
    hairShortData, hairShortData, hairShortData, hairShortData, hairShortData,
    hairShortData, hairShortData, hairShortData, hairShortData, hairShortData,
    hairShortData, hairShortData, hairShortData, hairShortData, hairShortData,
    hairShortData, hairShortData, hairShortData, hairShortData, hairShortData,
    hairShortData, hairShortData, hairShortData, hairShortData, hairShortData, hairShortData
);
hairClumpTransforms.push(
    { x: -0.21, y: 0.03, z: 0, rotX: 5, rotY: 100, scale: 0.4 },
    { x: -0.24, y: 0.03, z: 0.05, rotX: 2, rotY: 120, scale: 0.4 },
    { x: -0.24, y: 0.03, z: 0.1, rotX: 0, rotY: 130, scale: 0.7 },
    { x: -0.21, y: 0.03, z: 0.15, rotX: 0, rotY: 140, scale: 0.7 },
    { x: -0.18, y: 0.05, z: 0.18, rotX: 0, rotY: 150, scale: 0.7 },
    { x: -0.15, y: 0.05, z: 0.2, rotX: -10, rotY: 160, scale: 0.7 },
    { x: -0.07, y: 0.05, z: 0.2, rotX: -10, rotY: 170, scale: 0.55 },
    { x: -0.10, y: 0.05, z: 0.2, rotX: -10, rotY: 170, scale: 0.51 },
    { x: -0.03, y: 0.05, z: 0.2, rotX: -10, rotY: 170, scale: 0.52 },
    { x: 0.0, y: 0.08, z: 0.2, rotX: -10, rotY: 180, scale: 0.5 },
    { x: 0.03, y: 0.05, z: 0.2, rotX: -10, rotY: 190, scale: 0.52 },
    { x: 0.07, y: 0.05, z: 0.2, rotX: -10, rotY: 190, scale: 0.51 },
    { x: 0.10, y: 0.05, z: 0.2, rotX: -10, rotY: 190, scale: 0.55 },
    { x: 0.15, y: 0.05, z: 0.2, rotX: -10, rotY: 200, scale: 0.7 },
    { x: 0.18, y: 0.05, z: 0.18, rotX: 0, rotY: 210, scale: 0.7 },
    { x: 0.21, y: 0.03, z: 0.15, rotX: 0, rotY: 220, scale: 0.7 },
    { x: 0.24, y: 0.03, z: 0.1, rotX: 0, rotY: 230, scale: 0.7 },
    { x: 0.24, y: 0.03, z: 0.05, rotX: 2, rotY: 240, scale: 0.4 },
    { x: 0.21, y: 0.03, z: 0, rotX: 5, rotY: 260, scale: 0.4 },
    { x: 0.18, y: 0.02, z: -0.05, rotX: 10, rotY: 280, scale: 0.45 },
    { x: 0.12, y: 0.02, z: -0.08, rotX: 12, rotY: 300, scale: 0.5 },
    { x: 0.06, y: 0.02, z: -0.1, rotX: 15, rotY: 320, scale: 0.6 },
    { x: 0.00, y: 0.02, z: -0.12, rotX: 18, rotY: 340, scale: 0.6 },
    { x: -0.06, y: 0.02, z: -0.1, rotX: 15, rotY: 0, scale: 0.6 },
    { x: -0.12, y: 0.02, z: -0.08, rotX: 12, rotY: 20, scale: 0.5 },
    { x: -0.18, y: 0.02, z: -0.05, rotX: 10, rotY: 40, scale: 0.45 },
    { x: -0.20, y: 0.03, z: -0.03, rotX: 8, rotY: 60, scale: 0.4 },
    { x: -0.21, y: 0.03, z: -0.01, rotX: 5, rotY: 80, scale: 0.4 },
    { x: 0.21, y: 0.03, z: -0.01, rotX: 5, rotY: 260, scale: 0.4 },
    { x: 0.20, y: 0.03, z: -0.03, rotX: 8, rotY: 280, scale: 0.4 },
    { x: 0.18, y: 0.03, z: -0.05, rotX: 10, rotY: 300, scale: 0.45 },
    { x: 0.12, y: 0.03, z: -0.08, rotX: 12, rotY: 320, scale: 0.5 },
    { x: 0.06, y: 0.03, z: -0.1, rotX: 15, rotY: 340, scale: 0.55 }
);
hairClumpData.push(
    hairLongData, hairLongData, hairLongData, hairLongData, hairLongData,
    hairLongData, hairLongData, hairLongData, hairLongData
);
hairClumpTransforms.push(
    { x: -0.22, y: 0.02, z: 0.00, rotZ: -6, rotX: 12, rotY: 45, scale: 1.0 },
    { x: -0.17, y: 0.01, z: 0.00, rotZ: -5, rotX: 14, rotY: 30, scale: 1.0 },
    { x: -0.12, y: 0.05, z: 0.00, rotZ: -4, rotX: 16, rotY: 15, scale: 1.05 },
    { x: -0.06, y: 0.05, z: 0.00, rotZ: -3, rotX: 18, rotY: 5, scale: 1.05 },
    { x: 0.00, y: 0.05, z: 0.00, rotZ: 0, rotX: 20, rotY: 0, scale: 1.1 },
    { x: 0.06, y: 0.05, z: 0.00, rotZ: 3, rotX: 18, rotY: -5, scale: 1.05 },
    { x: 0.12, y: 0.05, z: 0.00, rotZ: 4, rotX: 16, rotY: -15, scale: 1.05 },
    { x: 0.17, y: 0.01, z: 0.00, rotZ: 5, rotX: 14, rotY: -30, scale: 1.0 },
    { x: 0.22, y: 0.02, z: 0.00, rotZ: 6, rotX: 12, rotY: -45, scale: 1.0 },
);
var allManualHairData = combineGeometry(hairClumpData, hairClumpTransforms);

var handFanData = combineGeometry(
    [leafData, leafDataMid, leafData],
    [
        { x: 0, y: 0, z: 0, rotZ: 40 },
        { x: 0, y: 0, z: 0, rotZ: 0 },
        { x: 0, y: 0, z: 0, rotZ: -40 }
    ]
);

var armDataLeft = combineGeometry(
    [armData, handFanData],
    [
        { x: 0.25, y: 0, z: 0, rotZ: -90 },
        { x: 0.3, y: 0, z: 0, rotZ: -90, rotY: -90 }
    ]
);

var armDataRight = combineGeometry(
    [armData, handFanData],
    [
        { x: -0.25, y: 0, z: 0, rotZ: 90 },
        { x: -0.3, y: 0, z: 0, rotZ: 90, rotY: 90 }
    ]
);

var oneLegData = combineGeometry(
    [thighData, legData, footData, heelData],
    [
        { x: 0, y: 0, z: 0 },
        { x: 0, y: -0.325, z: 0 },
        { x: 0, y: -0.45, z: 0 },
        { x: 0, y: -0.55, z: 0 }
    ]
);
var leftLegData = oneLegData;
var rightLegData = oneLegData;

var bodyData = combineGeometry(
    [
        torsoData, noseData,
        leftEyeData, rightEyeData,
        mouthData,
        toothData,
        toothData,
        outerEarData,
        innerEarData,
        outerEarData,
        innerEarData,
        hairBaseData,
        allManualHairData,
        leftLegData, rightLegData
    ],
    [
        { x: 0, y: 0, z: 0 }, // Torso
        { x: 0, y: 0.15, z: 0.45, rotX: 90 }, // Hidung
        { x: 0, y: 0, z: 0.04, rotY: -25 }, // Mata Kiri
        { x: 0, y: 0, z: 0.04, rotY: 25 },// Mata Kanan
        { x: 0, y: 0, z: 0.39 }, // Mulut
        { x: 0, y: 0.015, z: 0.4 }, // Gigi Atas
        { x: 0, y: -0.015, z: 0.4 }, // Gigi Bawah
        { x: -0.3, y: 1, z: 0.05, rotY: 90, rotZ: 15 }, // Telinga Luar Kiri
        { x: -0.3, y: 0.88, z: 0.06, rotY: 90, rotZ: 15 }, // Telinga Dalam Kiri
        { x: 0.3, y: 1, z: 0.05, rotY: 90, rotZ: -15 },  // Telinga Luar Kanan
        { x: 0.3, y: 0.88, z: 0.06, rotY: 90, rotZ: -15 }, // Telinga Dalam Kanan
        { x: 0, y: 0.25, z: 0.0 }, // Dasar Rambut
        { x: 0, y: 0.6, z: -0.05 }, // Rambut Manual
        { x: -0.4, y: -0.4, z: 0 }, // Kaki Kiri
        { x: 0.4, y: -0.4, z: 0 }, // Kaki Kanan
    ]
);


var razorLeafBlueprint = createDetailedLeaf(0.5, 0.15, 0.01, COLORS.LEAF_DARK, 4, 2);

var darkPulseBlueprint = createTorus(0.1, 0.04, 12, 8, COLORS.DARK_PULSE);

var leafStormBlueprint = razorLeafBlueprint;

var swordBladeData = createBox(0.1, 1.0, 0.02, COLORS.SWORD_BLUE);
var swordGuardData = createBox(0.3, 0.05, 0.05, COLORS.SWORD_HILT);
var swordHiltData = createBox(0.04, 0.2, 0.03, COLORS.SWORD_HILT);
var swordPommelData = createSphere(0.04, 6, 6, COLORS.SWORD_HILT);

var swordsDanceBlueprint = combineGeometry(
    [swordBladeData, swordGuardData, swordHiltData, swordPommelData],
    [
        { x: 0, y: 0.125, z: 0 },
        { x: 0, y: -0.4, z: 0 },
        { x: 0, y: -0.525, z: 0 },
        { x: 0, y: -0.64, z: 0 }
    ]
);

var activeLeaves = [];
var leafLifetime = 150;
var leafSpawnSpeed = 0.15;

var activePulses = [];
var pulseLifetime = 100;
var pulseSpeed = 0.05;

var activeStorm = [];
var stormLifetime = 200;
var stormSpawnRate = 3;
var stormTimer = 0;

var activeSwords = [];
var swordsDanceDuration = 180;
var defaultAmbient = [0.6, 0.6, 0.6];


window.shiftryGeometries = {
    body: bodyData,
    shoulder: shoulderData,
    armLeft: armDataLeft,
    armRight: armDataRight,
    razorLeaf: razorLeafBlueprint,
    darkPulse: darkPulseBlueprint,
    leafStorm: leafStormBlueprint,
    sword: swordsDanceBlueprint
};