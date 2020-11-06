#include"libOne.h"
//----------------------------------------------------------------------
// シェープの輪郭となる頂点位置を一筆書きとなるようにプログラムしていく
// 最初の頂点からトライアングルファンで塗りつぶせるように考慮する必要がある
//----------------------------------------------------------------------

//ひし形をつくる
int createDiamondShape() {
    //ひし形の対角線の半分の長さlx ly
    float lx = 0.5f;
    float ly = 0.7f;
    //頂点位置
    SHAPE_VERTEX vertices[] = {
        0, ly,
        lx, 0,
        0, -ly,
        -lx, 0,
    };
    //シェープを作って番号を返す
    return createShape(vertices, sizeof(vertices) / sizeof(SHAPE_VERTEX));
}

//星形をつくる
int createStarShape() {
    //頂点位置
    angleMode(DEGREES);
    const int numVertices = 10;
    SHAPE_VERTEX vertices[numVertices];
    float divDeg = 360.0f / numVertices;
    for (int i = 0; i < numVertices; i++) {
        float radius = 0.4f + 0.4f * (i % 2);
        float deg = divDeg * i;
        vertices[i].x = sin(deg) * radius;
        vertices[i].y = cos(deg) * radius;
    }
    //シェープを作って番号を返す
    return createShape(vertices, numVertices);
}


//角丸多角形
int createKadomaruShape(
    //扇型(角)の数
    int numCorners,
    //弧を分割する角度(輪郭の滑らかさ)
    int divDeg,
    //半径
    float radius,
    //扇型を中心から移動させる距離
    float vlen
) {
    angleMode(DEGREES);
    //扇型の中心角
    const int angle = 360 / numCorners;
    const int numVertices = (angle / divDeg + 1) * numCorners;
    //最初に頂点を用意する角度
    float offsetDeg = (180 - angle) / 2.0f;
    SHAPE_VERTEX* vertices = new SHAPE_VERTEX[numVertices];
    for (int i = 0; i < numVertices; i++) {
        int w = i / (numVertices / numCorners);
        float vx = -sin((float)angle * w) * vlen;
        float vy = -cos((float)angle * w) * vlen;
        float deg = offsetDeg + divDeg * (i - w);
        vertices[i].x = vx + cos(deg) * radius;
        vertices[i].y = vy + -sin(deg) * radius;
    }
    //　シェープを作って番号を返す
    int shapeIdx = createShape(vertices, numVertices);
    delete[] vertices;
    return shapeIdx;
}

//カプセルをつくる
int createCapsuleShape() {
    angleMode(DEGREES);
    //弧を分割する角度
    const int divDeg = 10;
    //半径
    float radius = 0.5f;
    const int numVertices = (180 / divDeg + 1) * 2;
    SHAPE_VERTEX vertices[numVertices];
    //上の半円
    float vy = -0.1f;
    for (int i = 0; i < numVertices / 2; i++) {
        float deg = (float)divDeg * i;
        vertices[i].x = cos(deg) * radius;
        vertices[i].y = vy + -sin(deg) * radius;
    }
    //下の半円
    vy = -vy;
    for (int i = numVertices / 2; i < numVertices; i++) {
        float deg = (float)divDeg * (i - 1);
        vertices[i].x = cos(deg) * radius;
        vertices[i].y = vy + -sin(deg) * radius;
    }

    //　シェープを作って番号を返す
    return createShape(vertices, numVertices);
}

//おにぎり（角丸正三角形）をつくる
int createOnigiriShape() {
    angleMode(DEGREES);
    //曲線を分割する角度
    const int divDeg = 10;
    const int numVertices = (120 / divDeg + 1) * 3;
    SHAPE_VERTEX vertices[numVertices];
    float offsetDeg = 30;
    float radius = 0.5f;
    float vl = 0.1f;
    for (int i = 0; i < numVertices; i++) {
        int w = i / (numVertices / 3);
        float vx = -sin(120.0f * w) * vl;
        float vy = -cos(120.0f * w) * vl;
        float deg = offsetDeg + divDeg * (i - w);
        vertices[i].x = vx + cos(deg) * radius;
        vertices[i].y = vy + -sin(deg) * radius;
    }
    /*
        //上の丸角
        float vx = -sin(0)*vl;
        float vy = -cos(0)*vl;
        for (int i = 0; i < numVertices / 3; i++) {
            float deg = offsetDeg + divDeg * i;
            vertices[i].x = vx + cos(deg) * radius;
            vertices[i].y = vy + -sin(deg) * radius;
        }
        //左下の丸角
        vx = -sin(120)*vl;
        vy = -cos(120)*vl;
        for (int i = numVertices / 3; i < numVertices / 3 * 2; i++) {
            float deg = offsetDeg + divDeg * (i - 1);
            vertices[i].x = vx + cos(deg) * radius;
            vertices[i].y = vy + -sin(deg) * radius;
        }
        //右下の丸角
        vx = -sin(240)*vl;
        vy = -cos(240)*vl;
        for (int i = numVertices / 3 * 2; i < numVertices; i++) {
            float deg = offsetDeg + divDeg * (i - 2);
            vertices[i].x = vx + cos(deg) * radius;
            vertices[i].y = vy + -sin(deg) * radius;
        }
    */
    //　シェープを作って番号を返す
    return createShape(vertices, numVertices);
}

//角丸正方形をつくる
int createKakumaruShape() {
    angleMode(DEGREES);
    //曲線を分割する角度
    const int divDeg = 10;
    //頂点数
    const int numVertices = (90 / divDeg + 1) * 4;
    SHAPE_VERTEX vertices[numVertices];
    //半径
    float radius = 0.3f;
    //中心から離す距離
    float vl = 0.3f;
    for (int i = 0; i < numVertices; i++) {
        int w = i / (numVertices / 4);
        float vx = cos(45.0f + 90 * w) * vl;
        float vy = -sin(45.0f + 90 * w) * vl;
        float deg = (float)divDeg * (i - w);
        vertices[i].x = vx + cos(deg) * radius;
        vertices[i].y = vy + -sin(deg) * radius;
    }
    /*
        //右上の丸角
        float vx = cos(45) * vl;
        float vy = -sin(45) * vl;
        for (int i = 0; i < numVertices / 4; i++) {
            float deg = (float)divDeg * i;
            vertices[i].x = vx + cos(deg) * radius;
            vertices[i].y = vy + -sin(deg) * radius;
        }
        //左上の丸角
        vx = cos(135) * vl;
        vy = -sin(135) * vl;
        for (int i = numVertices / 4; i < numVertices / 4 * 2; i++) {
            float deg = (float)divDeg * (i - 1);
            vertices[i].x = vx + cos(deg) * radius;
            vertices[i].y = vy + -sin(deg) * radius;
        }
        //左下の丸角
        vx = cos(225) * vl;
        vy = -sin(225) * vl;
        for (int i = numVertices / 4 * 2; i < numVertices / 4 * 3; i++) {
            float deg = (float)divDeg * (i - 2);
            vertices[i].x = vx + cos(deg) * radius;
            vertices[i].y = vy + -sin(deg) * radius;
        }
        //左下の丸角
        vx = cos(315) * vl;
        vy = -sin(315) * vl;
        for (int i = numVertices / 4 * 3; i < numVertices; i++) {
            int debug = i / (numVertices / 4);
            float deg = (float)divDeg * (i - 3);
            vertices[i].x = vx + cos(deg) * radius;
            vertices[i].y = vy + -sin(deg) * radius;
        }
    */
    //　シェープを作って番号を返す
    return createShape(vertices, numVertices);
}

//ハート型をつくる
int createHeartShape() {
    angleMode(DEGREES);
    const int divDeg = 5;
    const int numVertices = 360 / divDeg;
    SHAPE_VERTEX vertices[numVertices];
    float scale = 1.0f / 24;
    for (int i = 0; i < numVertices; i++) {
        float deg = (float)divDeg * i;
        vertices[i].x =
            16 * pow(sin(deg), 3) * scale;
        vertices[i].y =
            (-13 * cos(deg)
                + 5 * cos(2 * deg)
                + 2 * cos(3 * deg)
                + 1 * cos(4 * deg)) * scale;
    }
    //　シェープを作って番号を返す
    return createShape(vertices, numVertices);
}

void gmain() {
    window(1280, 720, full);
    const int n = 6;
    int shapes[n];
    COLOR colors[n];
    shapes[0] = createKadomaruShape(6, 10, 0.3f, 0.3f);
    //shapes[0] = createCapsuleShape();
    shapes[1] = createDiamondShape();
    shapes[2] = createOnigiriShape();
    shapes[3] = createStarShape();
    shapes[4] = createKakumaruShape();
    shapes[5] = createHeartShape();
    colors[0] = COLOR(255, 187, 187);
    colors[1] = COLOR(170, 221, 221);
    colors[2] = COLOR(153, 221, 255);
    colors[3] = COLOR(255, 255, 187);
    colors[4] = COLOR(221, 238, 170);
    colors[5] = COLOR(255, 187, 221);
    //描画パラメータ
    //縦横に並べる数
    float size = 15;
    stroke(128, 128, 128);
    float sw = 1;
    float deg = 0;
    bool rotateFlag = true;
    enum STATE { LEVEL1, LEVEL2, LEVEL3 };
    STATE state = STATE::LEVEL1;
    struct LEVEL {
        int nx, ny;
    };
    struct LEVEL level[3] = {
        6,1,
        6,3,
        8,5
    };
    int nx = level[state].nx;
    int ny = level[state].ny;
    while (notQuit) {
        switch (state) {
        case STATE::LEVEL1:
            if (isTrigger(KEY_D)) {
                state = STATE::LEVEL2;
                nx = level[state].nx;
                ny = level[state].ny;
            }
            break;
        case STATE::LEVEL2:
            if (isTrigger(KEY_A)) {
                state = STATE::LEVEL1;
                nx = level[state].nx;
                ny = level[state].ny;
            }
            if (isTrigger(KEY_D)) {
                state = STATE::LEVEL3;
                nx = level[state].nx;
                ny = level[state].ny;
            }
            break;
        case STATE::LEVEL3:
            if (isTrigger(KEY_D)) {
                nx *= 2;
                ny *= 2;
                if (nx > 128) {
                    nx = 128;
                    ny = 60;
                }
            }
            if (isTrigger(KEY_LEFT) || isTrigger(KEY_A)) {
                nx /= 2;
                ny /= 2;
                if (nx < 8) {
                    state = STATE::LEVEL2;
                    nx = level[state].nx;
                    ny = level[state].ny;
                }
            }
            break;
        }
        //シェイプ同士の間隔
        float dx = Width / nx;
        float dy = Height / ny;
        size = dx / 2;
        //輪郭線の太さ
        if (isPress(KEY_W)) {
            sw += 0.25f;
        }
        if (isPress(KEY_S)) {
            sw -= 0.25f;
            if (sw < 0.0f) { sw = 0.0f; }
        }
        strokeWeight(sw);
        //回転するしない
        if (isTrigger(KEY_SPACE)) {
            rotateFlag = !rotateFlag;
        }
        if (rotateFlag) {
            deg += 1;
        }
        else {
            deg = 0;
        }
        //描画
        clear(255, 221, 238);
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int idx = (j + i) % n;
                float x = dx / 2 + dx * i;
                float y = dy / 2 + dy * j;
                int dirRot = 1 - idx % 2 * 2;
                fill(colors[idx]);
                shape(shapes[idx], x, y, deg * dirRot, size);
            }
        }
    }
}