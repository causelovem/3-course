#include <iostream>
#include <vector>

#include "Utility.h"
#include "SOIL.h"

using namespace std;

const uint GRASS_INSTANCES = 10000; // Количество травинок

GL::Camera camera;               // Мы предоставляем Вам реализацию камеры. В OpenGL камера - это просто 2 матрицы. Модельно-видовая матрица и матрица проекции. // ###
                                 // Задача этого класса только в том чтобы обработать ввод с клавиатуры и правильно сформировать эти матрицы.
                                 // Вы можете просто пользоваться этим классом для расчёта указанных матриц.

GLuint grassPointsCount; // Количество вершин у модели травинки
GLuint grassShader;      // Шейдер, рисующий траву
GLuint grassVAO;         // VAO для травы (что такое VAO почитайте в доках)
GLuint grassVariance;    // Буфер для смещения координат травинок
vector<VM::vec4> grassVarianceData(GRASS_INSTANCES); // Вектор со смещениями для координат травинок
vector<VM::vec2> grassParamsData(GRASS_INSTANCES);

float wind_forceX = 0.05, wind_forceZ = 0.05, tau = 0.0, k = 5.0;

GLuint groundShader; // Шейдер для земли
GLuint groundVAO; // VAO для земли

int ground_x = 2, ground_z = 2;

GLuint ground_texture;
GLuint grass_texture;

// Размеры экрана
uint screenWidth = 800;
uint screenHeight = 600;

// Это для захвата мышки. Вам это не потребуется (это не значит, что нужно удалять эту строку)
bool captureMouse = true;

// Функция, рисующая замлю
void DrawGround()
{
    // Используем шейдер для земли
    glUseProgram(groundShader);                                                  CHECK_GL_ERRORS

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, ground_texture);
    glUniform1i(glGetUniformLocation(groundShader, "point"), 0);
    //glUniform1i(groundShader, 0);
    // Устанавливаем юниформ для шейдера. В данном случае передадим перспективную матрицу камеры
    // Находим локацию юниформа 'camera' в шейдере
    GLint cameraLocation = glGetUniformLocation(groundShader, "camera");         CHECK_GL_ERRORS
    // Устанавливаем юниформ (загружаем на GPU матрицу проекции?)                                                     // ###
    glUniformMatrix4fv(cameraLocation, 1, GL_TRUE, camera.getMatrix().data().data()); CHECK_GL_ERRORS

    // Подключаем VAO, который содержит буферы, необходимые для отрисовки земли
    glBindVertexArray(groundVAO);                                                CHECK_GL_ERRORS

    // Рисуем землю: 2 треугольника (6 вершин)
    glDrawArrays(GL_TRIANGLES, 0, 6);                                            CHECK_GL_ERRORS

    // Отсоединяем VAO
    glBindVertexArray(0);                                                        CHECK_GL_ERRORS
    // Отключаем шейдер
    glUseProgram(0);                                                             CHECK_GL_ERRORS
}

// Обновление смещения травинок
void UpdateGrassVariance()
{
    // Генерация случайных смещений
    /*X.. + kx = w_f*/
    tau += 0.01;

    for (uint i = 0; i < GRASS_INSTANCES / 2; ++i)
    {
        grassVarianceData[i].x = wind_forceX / k * (1 - cos(sqrt(k) * tau));
        grassVarianceData[i].z = wind_forceZ / k * (1 - cos(sqrt(k) * tau));
        //grassVarianceData[i].y = - wind_forceZ / k * (1 - cos(sqrt(k) * tau));
        grassVarianceData[i].y = - sqrt(grassVarianceData[i].x *  grassVarianceData[i].x + grassVarianceData[i].z * grassVarianceData[i].z);

        grassVarianceData[GRASS_INSTANCES / 2 + i].x = wind_forceX / k * (1 - cos(sqrt(k) * tau + 1));
        grassVarianceData[GRASS_INSTANCES / 2 + i].z = wind_forceZ / k * (1 - cos(sqrt(k) * tau + 1));
        //grassVarianceData[GRASS_INSTANCES / 2 + i].y = - wind_forceZ / k * (1 - cos(sqrt(k) * tau + 1));
        grassVarianceData[GRASS_INSTANCES / 2 + i].y = - sqrt(grassVarianceData[GRASS_INSTANCES / 2 + i].x * grassVarianceData[GRASS_INSTANCES / 2 + i].x + grassVarianceData[GRASS_INSTANCES / 2 + i].z * grassVarianceData[GRASS_INSTANCES / 2 + i].z);
    }

    // Привязываем буфер, содержащий смещения
    glBindBuffer(GL_ARRAY_BUFFER, grassVariance);                                CHECK_GL_ERRORS
    // Загружаем данные в видеопамять
    glBufferData(GL_ARRAY_BUFFER, sizeof(VM::vec4) * GRASS_INSTANCES, grassVarianceData.data(), GL_STATIC_DRAW); CHECK_GL_ERRORS
    // Отвязываем буфер
    glBindBuffer(GL_ARRAY_BUFFER, 0);                                            CHECK_GL_ERRORS
}

// Рисование травы
void DrawGrass()
{
    // Тут то же самое, что и в рисовании земли
    glUseProgram(grassShader);                                                   CHECK_GL_ERRORS

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, grass_texture);
    glUniform1i(glGetUniformLocation(grassShader, "point"), 0);

    GLint cameraLocation = glGetUniformLocation(grassShader, "camera");          CHECK_GL_ERRORS
    glUniformMatrix4fv(cameraLocation, 1, GL_TRUE, camera.getMatrix().data().data()); CHECK_GL_ERRORS

    glBindVertexArray(grassVAO);                                                 CHECK_GL_ERRORS
    // Обновляем смещения для травы
    UpdateGrassVariance();
    // Отрисовка травинок в количестве GRASS_INSTANCES
    glDrawArraysInstanced(GL_TRIANGLES, 0, grassPointsCount, GRASS_INSTANCES);   CHECK_GL_ERRORS
    glBindVertexArray(0);                                                        CHECK_GL_ERRORS
    glUseProgram(0);                                                             CHECK_GL_ERRORS
}

// Эта функция вызывается для обновления экрана
void RenderLayouts()
{
    // Включение буфера глубины
    glEnable(GL_DEPTH_TEST);
    //glClearColor(135.0 / 255.0, 206.0 / 255.0, 235.0 / 255.0, 0.0);
    glClearColor(0.0 / 255.0, 191.0 / 255.0, 255.0 / 255.0, 0.0);
    // Очистка буфера глубины и цветового буфера
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Рисуем меши
    DrawGround();
    DrawGrass();  
    glutSwapBuffers();
}

// Завершение программы
void FinishProgram()
{
    glutDestroyWindow(glutGetWindow());
}

// Обработка события нажатия клавиши (специальные клавиши обрабатываются в функции SpecialButtons)
void KeyboardEvents(unsigned char key, int x, int y)
{
    if (key == 27) {
        FinishProgram();
    } else if (key == 'w') {
        camera.goForward();
    } else if (key == 's') {
        camera.goBack();
    } else if (key == 'm') {
        captureMouse = !captureMouse;
        if (captureMouse) {
            glutWarpPointer(screenWidth / 2, screenHeight / 2);
            glutSetCursor(GLUT_CURSOR_NONE);
        } else {
            glutSetCursor(GLUT_CURSOR_RIGHT_ARROW);
        }
    }
}

// Обработка события нажатия специальных клавиш
void SpecialButtons(int key, int x, int y)
{
    if (key == GLUT_KEY_RIGHT) {
        camera.rotateY(0.02);
    } else if (key == GLUT_KEY_LEFT) {
        camera.rotateY(-0.02);
    } else if (key == GLUT_KEY_UP) {
        camera.rotateTop(-0.02);
    } else if (key == GLUT_KEY_DOWN) {
        camera.rotateTop(0.02);
    }
}

void IdleFunc()
{
    glutPostRedisplay();
}

// Обработка события движения мыши
void MouseMove(int x, int y)
{
    if (captureMouse) {
        int centerX = screenWidth / 2,
            centerY = screenHeight / 2;
        if (x != centerX || y != centerY) {
            camera.rotateY((x - centerX) / 1000.0f);
            camera.rotateTop((y - centerY) / 1000.0f);
            glutWarpPointer(centerX, centerY);
        }
    }
}

// Обработка нажатия кнопки мыши
void MouseClick(int button, int state, int x, int y)
{
}

// Событие изменение размера окна
void windowReshapeFunc(GLint newWidth, GLint newHeight)
{
    glViewport(0, 0, newWidth, newHeight);
    screenWidth = newWidth;
    screenHeight = newHeight;

    camera.screenRatio = (float)screenWidth / screenHeight;
}

// Инициализация окна
void InitializeGLUT(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE); // GLUT_DEPTH
    glutInitContextVersion(3, 0);
    glutInitWindowPosition(-1, -1);
    glutInitWindowSize(screenWidth, screenHeight);
    glutCreateWindow("Computer Graphics 3");
    glutWarpPointer(400, 300);
    glutSetCursor(GLUT_CURSOR_NONE);

    glutDisplayFunc(RenderLayouts);
    glutKeyboardFunc(KeyboardEvents);
    glutSpecialFunc(SpecialButtons);
    glutIdleFunc(IdleFunc);
    glutPassiveMotionFunc(MouseMove);
    glutMouseFunc(MouseClick);
    glutReshapeFunc(windowReshapeFunc);
}

// Генерация позиций травинок (эту функцию вам придётся переписать)
vector<VM::vec2> GenerateGrassPositions()
{
    vector<VM::vec2> grassPositions(GRASS_INSTANCES);
    for (uint i = 0; i < GRASS_INSTANCES; ++i)
    {
        grassPositions[i] = VM::vec2((float)rand() * ground_x / RAND_MAX, (float)rand() * ground_z / RAND_MAX);
        grassParamsData[i] = VM::vec2((float)i, float(float(rand() % 8 + 5) / 100.0));
    }

    return grassPositions;
}

// Здесь вам нужно будет генерировать меш 
vector<VM::vec4> GenMesh(uint n)
{
    /*return {
        // 1 tringle
        VM::vec4(0, 0, 0, 1),
        VM::vec4(1, 0, 0, 1),
        VM::vec4(1, 1.0 / 3.0, 0, 1),

        // 2 tringle
        VM::vec4(0, 0, 0, 1),
        VM::vec4(0, 1.0 / 3.0, 0, 1),
        VM::vec4(1, 1.0 / 3.0, 0, 1),

        // 3 tringle
        VM::vec4(0, 1.0 / 3.0, 0, 1),
        VM::vec4(1, 1.0 / 3.0, 0, 1),
        VM::vec4(17.0 / 20.0, 2.0 / 3.0, 0, 1),

        // 4 tringle
        VM::vec4(0, 1.0 / 3.0, 0, 1),
        VM::vec4(3.0 / 20.0, 2.0 / 3.0, 0, 1),
        VM::vec4(17.0 / 20.0, 2.0 / 3.0, 0, 1),

        // 5 tringle
        VM::vec4(3.0 / 20.0, 2.0 / 3.0, 0, 1),
        VM::vec4(17.0 / 20.0, 2.0 / 3.0, 0, 1),
        VM::vec4(1.0 / 2.0, 1, 0, 1),
    };*/

    return {
        // 1 tringle
        VM::vec4(0, 0, 0, 1),
        VM::vec4(1, 0, 0, 1),
        VM::vec4(1, 1.0 / 4.0, 0, 1),

        // 2 tringle
        VM::vec4(0, 0, 0, 1),
        VM::vec4(0, 1.0 / 4.0, 0, 1),
        VM::vec4(1, 1.0 / 4.0, 0, 1),

        // 3 tringle
        VM::vec4(0, 1.0 / 4.0, 0, 1),
        VM::vec4(1, 1.0 / 4.0, 0, 1),
        VM::vec4(90.0 / 100.0, 2.0 / 4.0, 0, 1),

        // 4 tringle 
        VM::vec4(0, 1.0 / 4.0, 0, 1),
        VM::vec4(10.0 / 100.0, 2.0 / 4.0, 0, 1),
        VM::vec4(90.0 / 100.0, 2.0 / 4.0, 0, 1),

        // 5 tringle
        VM::vec4(10.0 / 100.0, 2.0 / 4.0, 0, 1),
        VM::vec4(90.0 / 100.0, 2.0 / 4.0, 0, 1),
        VM::vec4(80.0 / 100.0, 3.0 / 4.0, 0, 1),

        // 6 tringle
        VM::vec4(10.0 / 100.0, 2.0 / 4.0, 0, 1),
        VM::vec4(20.0 / 100.0, 3.0 / 4.0, 0, 1),
        VM::vec4(80.0 / 100.0, 3.0 / 4.0, 0, 1),

        // 7 tringle
        VM::vec4(20.0 / 100.0, 3.0 / 4.0, 0, 1),
        VM::vec4(80.0 / 100.0, 3.0 / 4.0, 0, 1),
        VM::vec4(50.0 / 100.0, 1, 0, 1),
    };
}

// Создание травы
void CreateGrass()
{
    //Текстура
    grass_texture = SOIL_load_OGL_texture("../Texture/grass.jpg",
    SOIL_LOAD_AUTO,
    SOIL_CREATE_NEW_ID,
    SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT
    );

    glBindTexture(GL_TEXTURE_2D, grass_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);                CHECK_GL_ERRORS
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);                CHECK_GL_ERRORS
    // Set texture filtering
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);            CHECK_GL_ERRORS
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);            CHECK_GL_ERRORS

    uint LOD = 1;
    // Создаём меш
    vector<VM::vec4> grassPoints = GenMesh(LOD);
    // Сохраняем количество вершин в меше травы
    grassPointsCount = grassPoints.size();
    // Создаём позиции для травинок
    vector<VM::vec2> grassPositions = GenerateGrassPositions();
    // Инициализация смещений для травинок
    for (uint i = 0; i < GRASS_INSTANCES; ++i) {
        grassVarianceData[i] = VM::vec4(0, 0, 0, 0);
    }

    /* Компилируем шейдеры
    Эта функция принимает на вход название шейдера 'shaderName',
    читает файлы shaders/{shaderName}.vert - вершинный шейдер
    и shaders/{shaderName}.frag - фрагментный шейдер,
    компилирует их и линкует.
    */
    grassShader = GL::CompileShaderProgram("grass");

    // Здесь создаём буфер
    GLuint pointsBuffer;
    // Это генерация одного буфера (в pointsBuffer хранится идентификатор буфера)
    glGenBuffers(1, &pointsBuffer);                                              CHECK_GL_ERRORS
    // Привязываем сгенерированный буфер
    glBindBuffer(GL_ARRAY_BUFFER, pointsBuffer);                                 CHECK_GL_ERRORS
    // Заполняем буфер данными из вектора
    glBufferData(GL_ARRAY_BUFFER, sizeof(VM::vec4) * grassPoints.size(), grassPoints.data(), GL_STATIC_DRAW); CHECK_GL_ERRORS

    // Создание VAO
    // Генерация VAO
    glGenVertexArrays(1, &grassVAO);                                             CHECK_GL_ERRORS
    // Привязка VAO
    glBindVertexArray(grassVAO);                                                 CHECK_GL_ERRORS

    // Получение локации параметра 'point' в шейдере
    GLuint pointsLocation = glGetAttribLocation(grassShader, "point");           CHECK_GL_ERRORS
    // Подключаем массив атрибутов к данной локации
    glEnableVertexAttribArray(pointsLocation);                                   CHECK_GL_ERRORS
    // Устанавливаем параметры для получения данных из массива (по 4 значение типа float на одну вершину)
    glVertexAttribPointer(pointsLocation, 4, GL_FLOAT, GL_FALSE, 0, 0);          CHECK_GL_ERRORS

    // Создаём буфер для позиций травинок
    GLuint positionBuffer;
    glGenBuffers(1, &positionBuffer);                                            CHECK_GL_ERRORS
    // Здесь мы привязываем новый буфер, так что дальше вся работа будет с ним до следующего вызова glBindBuffer
    glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);                               CHECK_GL_ERRORS
    glBufferData(GL_ARRAY_BUFFER, sizeof(VM::vec2) * grassPositions.size(), grassPositions.data(), GL_STATIC_DRAW); CHECK_GL_ERRORS

    GLuint positionLocation = glGetAttribLocation(grassShader, "position");      CHECK_GL_ERRORS
    glEnableVertexAttribArray(positionLocation);                                 CHECK_GL_ERRORS
    glVertexAttribPointer(positionLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);        CHECK_GL_ERRORS
    // Здесь мы указываем, что нужно брать новое значение из этого буфера для каждого инстанса (для каждой травинки)
    glVertexAttribDivisor(positionLocation, 1);                                  CHECK_GL_ERRORS

    // Создаём буфер для смещения травинок
    glGenBuffers(1, &grassVariance);                                             CHECK_GL_ERRORS
    glBindBuffer(GL_ARRAY_BUFFER, grassVariance);                                CHECK_GL_ERRORS
    glBufferData(GL_ARRAY_BUFFER, sizeof(VM::vec4) * GRASS_INSTANCES, grassVarianceData.data(), GL_STATIC_DRAW); CHECK_GL_ERRORS

    GLuint varianceLocation = glGetAttribLocation(grassShader, "variance");      CHECK_GL_ERRORS
    glEnableVertexAttribArray(varianceLocation);                                 CHECK_GL_ERRORS
    glVertexAttribPointer(varianceLocation, 4, GL_FLOAT, GL_FALSE, 0, 0);        CHECK_GL_ERRORS
    glVertexAttribDivisor(varianceLocation, 1);                                  CHECK_GL_ERRORS

    // Создаём буфер для поворота и высоты травинок
    GLuint paramsBuffer;
    glGenBuffers(1, &paramsBuffer);                                              CHECK_GL_ERRORS
    glBindBuffer(GL_ARRAY_BUFFER, paramsBuffer);                                 CHECK_GL_ERRORS
    glBufferData(GL_ARRAY_BUFFER, sizeof(VM::vec2) * GRASS_INSTANCES, grassParamsData.data(), GL_STATIC_DRAW); CHECK_GL_ERRORS

    GLuint paramsLocation = glGetAttribLocation(grassShader, "params");          CHECK_GL_ERRORS
    glEnableVertexAttribArray(paramsLocation);                                   CHECK_GL_ERRORS
    glVertexAttribPointer(paramsLocation, 2, GL_FLOAT, GL_FALSE, 0, 0);          CHECK_GL_ERRORS
    glVertexAttribDivisor(paramsLocation, 1);                                    CHECK_GL_ERRORS

    // Отвязываем VAO
    glBindVertexArray(0);                                                        CHECK_GL_ERRORS
    // Отвязываем буфер
    glBindBuffer(GL_ARRAY_BUFFER, 0);                                            CHECK_GL_ERRORS
}

// Создаём камеру (Если шаблонная камера вам не нравится, то можете переделать, но я бы не стал)
void CreateCamera()
{
    camera.angle = 45.0f / 180.0f * M_PI;
    camera.direction = VM::vec3(0, 0.3, -1);
    camera.position = VM::vec3(0.5, 0.2, 0);
    camera.screenRatio = (float)screenWidth / screenHeight;
    camera.up = VM::vec3(0, 1, 0);
    camera.zfar = 50.0f;
    camera.znear = 0.05f;
}

// Создаём замлю
void CreateGround()
{
    //Текстура
    ground_texture = SOIL_load_OGL_texture("../Texture/ground.jpg",
    SOIL_LOAD_AUTO,
    SOIL_CREATE_NEW_ID,
    SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT
    );

    glBindTexture(GL_TEXTURE_2D, ground_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);                CHECK_GL_ERRORS
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);                CHECK_GL_ERRORS
    // Set texture filtering
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);            CHECK_GL_ERRORS
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);            CHECK_GL_ERRORS
    
    // Земля состоит из двух треугольников
    vector<VM::vec4> meshPoints = {
        VM::vec4(0, 0, 0, 1),
        VM::vec4(ground_x, 0, 0, 1),
        VM::vec4(ground_x, 0, ground_z, 1),
        VM::vec4(0, 0, 0, 1),
        VM::vec4(ground_x, 0, ground_z, 1),
        VM::vec4(0, 0, ground_z, 1),
    };

    // Подробнее о том, как это работает читайте в функции CreateGrass

    groundShader = GL::CompileShaderProgram("ground");

    GLuint pointsBuffer;
    glGenBuffers(1, &pointsBuffer);                                              CHECK_GL_ERRORS
    glBindBuffer(GL_ARRAY_BUFFER, pointsBuffer);                                 CHECK_GL_ERRORS
    glBufferData(GL_ARRAY_BUFFER, sizeof(VM::vec4) * meshPoints.size(), meshPoints.data(), GL_STATIC_DRAW); CHECK_GL_ERRORS

    glGenVertexArrays(1, &groundVAO);                                            CHECK_GL_ERRORS
    glBindVertexArray(groundVAO);                                                CHECK_GL_ERRORS

    GLuint index = glGetAttribLocation(groundShader, "point");                   CHECK_GL_ERRORS
    glEnableVertexAttribArray(index);                                            CHECK_GL_ERRORS
    glVertexAttribPointer(index, 4, GL_FLOAT, GL_FALSE, 0, 0);                   CHECK_GL_ERRORS

    glBindVertexArray(0);                                                        CHECK_GL_ERRORS
    glBindBuffer(GL_ARRAY_BUFFER, 0);                                            CHECK_GL_ERRORS
}

int main(int argc, char **argv)
{
    srand(time(0));

    //wind_force = (uint)rand() % 5;

    putenv("MESA_GL_VERSION_OVERRIDE=3.3COMPAT");
    try {
        cout << "Start" << endl;
        InitializeGLUT(argc, argv);
        cout << "GLUT inited" << endl;
        glewInit();
        cout << "glew inited" << endl;
        CreateCamera();
        cout << "Camera created" << endl;
        CreateGrass();
        cout << "Grass created" << endl;
        CreateGround();
        cout << "Ground created" << endl;
        glutMainLoop();
    } catch (string s) {
        cout << s << endl;
    }
}
