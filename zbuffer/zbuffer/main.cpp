// Std. Includes
#include <string>
#include <algorithm>
using namespace std;

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// GLFW
#include <GLFW/glfw3.h>

// GL includes
#include "learnopengl/Shader.h"
#include "learnopengl/Camera.h"

// GLM Mathemtics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

//assimp
#include <assimp/Importer.hpp>      // 导入器在该头文件中定义
#include <assimp/scene.h>           // 读取到的模型数据都放在scene中
#include <assimp/postprocess.h>     // 该头文件中包含后处理的标志位定义

// Other Libs
#include "includes\FreeImage\TextureManager.h"
#include "Triple.h"
// Properties
GLuint screenWidth = 800, screenHeight = 600;
vector<Triple<GLfloat>> vertex;
vector<vector<Triple<GLfloat>>> face;
// Function prototypes

// Camera
Camera camera(glm::vec3(0.0f, 0.0f, 3.0f));
bool keys[1024];
GLfloat lastX = 400, lastY = 300;
bool firstMouse = true;
GLfloat scale=50.0f;
GLfloat *pixels;//帧缓存

GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;
GLint g_polygon_id=0;
struct Point {
	GLint x, y;
	GLfloat z;
};

struct ClassifiedPolygon{
	GLdouble a, b, c, d;
	GLint id;
	GLint dy;
	Triple<GLfloat> color;
	ClassifiedPolygon* next;
};

struct  ClassifiedEdge{
	GLdouble x;
	GLdouble dx;
	GLint dy;
	GLint id;
	bool used;
	ClassifiedEdge* next;
};

struct ActivePolygon {
	GLdouble a, b, c, d;
	GLint id;
	GLint dy;
	Triple<GLfloat> color;
	ActivePolygon* next;
};

struct ActiveEdgePair{
	GLdouble xl;
	GLdouble dxl;
	GLint dyl;
	GLdouble xr;
	GLdouble dxr;
	GLint dyr;
	GLdouble zl;
	GLdouble dzx;
	GLdouble dzy;
	GLint id;
	Triple<GLfloat> color;
	ActiveEdgePair* next;
};

vector<ClassifiedPolygon*> tClassifiedPolygon;

vector<ClassifiedEdge*> tClassifiedEdge;

ActivePolygon tActivePolygonHead;

ActiveEdgePair tActiveEdgePairHead;
GLfloat g_zbuffer[2160];


void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void Do_Movement();
void zbuffer();
void setBGColor(GLint y);
void setPixColor(GLint x, GLint y, Triple<GLfloat> color);
void init();
void findEdge(ClassifiedEdge** e, GLint y, GLint polygon_id);
void activateNewPolygon(GLint y);
void deepthUpdate(GLint y);
void activeEdgeTableUpdate(GLint y);
void activePolygonTableUpdate();
void clearData();
void processNode(aiNode* node, const aiScene* scene);
void loadModel(string path);
Point roundVertex(Triple<GLfloat> v);
Triple<GLfloat> getPolygonColor(vector<GLdouble> coffs);
GLuint generateAttachmentTexture(GLboolean depth, GLboolean stencil);
GLfloat cubeVertices[] = {
	// Positions          // Texture Coords
	-0.5f, -0.5f, -0.5f,  0.0f, 0.0f,
	0.5f, -0.5f, -0.5f,  1.0f, 0.0f,
	0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
	0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
	-0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
	-0.5f, -0.5f, -0.5f,  0.0f, 0.0f,

	-0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
	0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
	0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
	0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
	-0.5f,  0.5f,  0.5f,  0.0f, 1.0f,
	-0.5f, -0.5f,  0.5f,  0.0f, 0.0f,

	-0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
	-0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
	-0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
	-0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
	-0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
	-0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

	0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
	0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
	0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
	0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
	0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
	0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

	-0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
	0.5f, -0.5f, -0.5f,  1.0f, 1.0f,
	0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
	0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
	-0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
	-0.5f, -0.5f, -0.5f,  0.0f, 1.0f,

	-0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
	0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
	0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
	0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
	-0.5f,  0.5f,  0.5f,  0.0f, 0.0f,
	-0.5f,  0.5f, -0.5f,  0.0f, 1.0f
};
// The MAIN function, from here we start our application and run our Game loop
int main()
{
	// Init GLFW
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

	GLFWwindow* window = glfwCreateWindow(screenWidth, screenHeight, "LearnOpenGL", nullptr, nullptr); // Windowed
	glfwMakeContextCurrent(window);

	// Set the required callback functions
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, mouse_callback);
	glfwSetScrollCallback(window, scroll_callback);

	// Options
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	// Initialize GLEW to setup the OpenGL Function pointers
	glewExperimental = GL_TRUE;
	glewInit();

	// Define the viewport dimensions
	glViewport(0, 0, screenWidth, screenHeight);


	// Setup and compile our shaders
	Shader screenShader("framebuffers_screen.vs", "framebuffers_screen.frag");

#pragma region "object_initialization"
	// Set the object data (buffers, vertex attributes)
	

	GLfloat quadVertices[] = {   // Vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
		// Positions   // TexCoords
		-1.0f,  1.0f,  0.0f, 1.0f,
		-1.0f, -1.0f,  0.0f, 0.0f,
		1.0f, -1.0f,  1.0f, 0.0f,

		-1.0f,  1.0f,  0.0f, 1.0f,
		1.0f, -1.0f,  1.0f, 0.0f,
		1.0f,  1.0f,  1.0f, 1.0f
	};

	// Setup screen VAO
	GLuint quadVAO, quadVBO;
	glGenVertexArrays(1, &quadVAO);
	glGenBuffers(1, &quadVBO);
	glBindVertexArray(quadVAO);
	glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (GLvoid*)0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (GLvoid*)(2 * sizeof(GLfloat)));
	glBindVertexArray(0);

#pragma endregion

	// Framebuffers
	GLuint framebuffer;
	glGenFramebuffers(1, &framebuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	pixels = (GLfloat*)malloc(screenHeight *screenWidth *3*sizeof(GLfloat));
	for (GLint i = 0; i < screenHeight * screenWidth;i++) {
		pixels[i * 3] = 0.2f;
		pixels[i * 3 + 1] = 1.0f;
		pixels[i * 3 + 2] = 0.2f;
	}
	GLuint pix_tex;
	glGenTextures(1, &pix_tex);
	// Draw as wireframe
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	// Game loop
	while (!glfwWindowShouldClose(window))
	{

		// Set frame time
		GLfloat currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;
		cout << 1/ deltaTime << endl;
		//glEnable(GL_DEPTH_TEST);

		/////////////////////////////////////////////////////
		// Bind to framebuffer and draw to color texture 
		// as we normally would.
		// //////////////////////////////////////////////////
		glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
		// Clear all attached buffers        
		glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // We're not using stencil buffer so why bother with clearing?

		/////////////////////////////////////////////////////
		// Bind to default framebuffer again and draw the 
		// quad plane with attched screen texture.
		// //////////////////////////////////////////////////
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		// Clear all relevant buffers
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set clear color to white (not really necessery actually, since we won't be able to see behind the quad anyways)
		glClear(GL_COLOR_BUFFER_BIT);
		glDisable(GL_DEPTH_TEST); // We don't care about depth information when rendering a single quad

								  // Draw Screen
		screenShader.Use();
		glBindVertexArray(quadVAO);
		init();//建立分类多边形表和分类边表
		zbuffer();
		clearData();
		glBindTexture(GL_TEXTURE_2D, pix_tex);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 800, 600, 0, GL_RGB, GL_FLOAT, pixels);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		// Set texture filtering
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		// Check and call events
		glfwPollEvents();
		Do_Movement();
		//glBindTexture(GL_TEXTURE_2D, textureColorbuffer);	
		// Use the color attachment texture as the texture of the quad plane
		//TextureManager::Inst()->BindTexture(cubeTexture);
		glBindTexture(GL_TEXTURE_2D,pix_tex);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		glBindVertexArray(0);

		// Swap the buffers
		glfwSwapBuffers(window);

	}

	// Clean up
	glDeleteFramebuffers(1, &framebuffer);
	glfwTerminate();
	return 0;
}

// This function loads a texture from file. Note: texture loading functions like these are usually 
// managed by a 'Resource Manager' that manages all resources (like textures, models, audio). 
// For learning purposes we'll just define it as a utility function.



// Generates a texture that is suited for attachments to a framebuffer
GLuint generateAttachmentTexture(GLboolean depth, GLboolean stencil)
{
	// What enum to use?
	GLenum attachment_type;
	if (!depth && !stencil)
		attachment_type = GL_RGB;
	else if (depth && !stencil)
		attachment_type = GL_DEPTH_COMPONENT;
	else if (!depth && stencil)
		attachment_type = GL_STENCIL_INDEX;

	//Generate texture ID and load texture data 
	GLuint textureID;
	glGenTextures(1, &textureID);
	glBindTexture(GL_TEXTURE_2D, textureID);
	if (!depth && !stencil)
		glTexImage2D(GL_TEXTURE_2D, 0, attachment_type, screenWidth, screenHeight, 0, attachment_type, GL_UNSIGNED_BYTE, NULL);
	else // Using both a stencil and depth test, needs special format arguments
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, screenWidth, screenHeight, 0, GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);

	return textureID;
}

#pragma region "User input"

// Moves/alters the camera positions based on user input
void Do_Movement()
{
	// Camera controls
	if (keys[GLFW_KEY_W])
		camera.ProcessKeyboard(FORWARD, deltaTime);
	if (keys[GLFW_KEY_S])
		camera.ProcessKeyboard(BACKWARD, deltaTime);
	if (keys[GLFW_KEY_A])
		camera.ProcessKeyboard(LEFT, deltaTime);
	if (keys[GLFW_KEY_D])
		camera.ProcessKeyboard(RIGHT, deltaTime);
}

// Is called whenever a key is pressed/released via GLFW
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	if (key==GLFW_KEY_Q) {
		scale += 2;
	}
	if (key==GLFW_KEY_E) {
		scale -= 2;
	}

	if (action == GLFW_PRESS)
		keys[key] = true;
	else if (action == GLFW_RELEASE)
		keys[key] = false;
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
	if (firstMouse)
	{
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	GLfloat xoffset = xpos - lastX;
	GLfloat yoffset = lastY - ypos;

	lastX = xpos;
	lastY = ypos;

	camera.ProcessMouseMovement(xoffset, yoffset);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	camera.ProcessMouseScroll(yoffset);
}

void zbuffer() {
	tActivePolygonHead.next = NULL;
	tActiveEdgePairHead.next = NULL;
	for (GLint y = screenHeight - 1; y >= 0;y--) {
		for (int i = 0; i < screenWidth;i++) {
			g_zbuffer[i] = -999999.0f;
		}//深度缓存重置为负无限大
		setBGColor(y);//重置帧缓存为背景色
		activateNewPolygon(y);//添加新的多边形到活化多边形表  
		deepthUpdate(y);//增量式深度更新帧缓存  
		activeEdgeTableUpdate(y);//活化边表元素修改  
		activePolygonTableUpdate();//更新活化多边形表  
	}
}

void setBGColor(GLint y) {
	for (int j = 0; j < screenWidth; j++) {
		pixels[(y * screenWidth + j) * 3] = 0.2f;
		pixels[(y * screenWidth + j) * 3 + 1] = 0.2f;
		pixels[(y * screenWidth + j) * 3 + 2] = 0.2f;
	}
}
void setPixColor(GLint x,GLint y,Triple<GLfloat> color) {
	pixels[(y * screenWidth + x) * 3] = color.x;
	pixels[(y * screenWidth + x) * 3 + 1] = color.y;
	pixels[(y * screenWidth + x) * 3 + 2] = color.z;
}


void init() {
	glm::mat4 model;
	glm::mat4 view = camera.GetViewMatrix();
	glm::mat4 projection = glm::perspective(camera.Zoom, (float)screenWidth / (float)screenHeight, 0.1f, 100.0f);
	model = glm::translate(model, glm::vec3(300.0f, 400.0f, 0.0f));
	GLfloat angle = (GLfloat)glfwGetTime();
	model = glm::rotate(model, angle, glm::vec3(1.0f, 1.0f, 0.0f));
	glm::mat4 MVP = view * model;
	for (int i = 0; i < screenWidth; i++) {
		g_zbuffer[i] = -999999.0f;
	}

	tClassifiedPolygon.clear();
	tClassifiedPolygon.assign(screenHeight, NULL);

	tClassifiedEdge.clear();
	tClassifiedEdge.assign(screenHeight, NULL);
	vertex.clear();
	face.clear();
	loadModel("model/planet.obj");
	g_polygon_id = 0;
	for (GLint i = 0; i < face.size();i++) {
		g_polygon_id++;
		vector<GLdouble> coefficient(4);
		coefficient[0] = (face[i][2].y - face[i][0].y)*(face[i][2].z - face[i][0].z) - (face[i][1].z - face[i][0].z)*(face[i][2].y - face[i][0].y);
		coefficient[1] = (face[i][2].x - face[i][0].x)*(face[i][1].z - face[i][0].z) - (face[i][1].x - face[i][0].x)*(face[i][2].z - face[i][0].z);
		coefficient[2] = (face[i][1].x - face[i][0].x)*(face[i][2].y - face[i][0].y) - (face[i][2].x - face[i][0].x)*(face[i][1].y - face[i][0].y);
		coefficient[3] = 0 - (coefficient[0] * face[i][0].x + coefficient[1] * face[i][0].y + coefficient[2] * face[i][0].z);
		//计算出a，b，c，d
		if (coefficient[2]<1e-8&&coefficient[2]>-1e-8)	continue;//垂直于xoy面的不考虑
		int polygon_maxy = (face[i][0].y > 0 ? (face[i][0].y + ((face[i][0].y < 0) ? -0.5 : 0.5)): 0), polygon_miny = polygon_maxy;
		for (GLint j = 0; j < face[i].size();j++) {
			Point a = roundVertex( face[i][j]);
			Point b = roundVertex(face[i][(j + 1) % face[i].size()]);
			if (a.y<b.y) {
				Point temp=a;
				a = b;
				b = temp;
			}//让a作为上顶点
			if (a.y < 0) continue;
			if (a.y - b.y < 1e-8 && a.y - b.y > -1e-8) continue;
			ClassifiedEdge* CE = (ClassifiedEdge*)malloc(sizeof(ClassifiedEdge));
			CE->x = a.x;
			CE->dx = (b.x - a.x + .0) / (a.y - b.y);
			CE->dy = a.y - b.y + 1;
			CE->id = g_polygon_id;
			CE->used = false;
			CE->next = NULL;
			if (a.y>=tClassifiedEdge.size())continue;
			if (tClassifiedEdge[a.y]==NULL) {
				tClassifiedEdge[a.y] = CE;
			}
			else{
				ClassifiedEdge* p_tmp = tClassifiedEdge[a.y];
				while (p_tmp->next)	p_tmp = p_tmp->next;
				p_tmp->next = CE;
			}
			if (a.y > polygon_maxy) polygon_maxy = a.y;
			if (b.y < polygon_miny) polygon_miny = b.y;
		}

		ClassifiedPolygon* CP = (ClassifiedPolygon*)malloc(sizeof(ClassifiedPolygon));
		CP->a = coefficient[0];
		CP->b = coefficient[1];
		CP->c = coefficient[2];
		CP->d = coefficient[3];
		CP->dy = polygon_maxy - polygon_miny + 1;
		CP->id = g_polygon_id;
		CP->color = getPolygonColor(coefficient);
		CP->next = NULL;
		if (polygon_maxy >= tClassifiedPolygon.size())	continue;
		if (tClassifiedPolygon[polygon_maxy]==NULL) {
			tClassifiedPolygon[polygon_maxy] = CP;
		}
		else{
			ClassifiedPolygon *p_tmp = tClassifiedPolygon[polygon_maxy];
			while (p_tmp->next)	p_tmp = p_tmp->next;
			p_tmp->next = CP;
		}
	}
}

Triple<GLfloat> getPolygonColor(vector<GLdouble> coffs) {
	//根据平面法向量与z轴（不分正反向）的夹角
	//夹角由小到大，颜色由明到暗
	if (coffs[0] < 0) coffs[0] = 0 - coffs[0];
	if (coffs[1] < 0) coffs[1] = 0 - coffs[1];
	if (coffs[2] < 0) coffs[2] = 0 - coffs[2];
	GLdouble costheta = coffs[2] / sqrt(coffs[0] + coffs[1] + coffs[2]);
	return Triple<GLfloat>(255 * costheta, 255 * costheta, 255 * costheta);
}

void activateNewPolygon(GLint y) {
	if (y < 0 || y >= tClassifiedPolygon.size())	return;
	if (tClassifiedPolygon[y]) {//将当前扫描线相关的多边形加入到活化多边形表
		ClassifiedPolygon* CP = tClassifiedPolygon[y];
		while (CP){
			ClassifiedPolygon *CP_tmp = (ClassifiedPolygon*)malloc(sizeof(ClassifiedPolygon));
			*CP_tmp = *CP;
			CP = CP_tmp;
			CP_tmp = CP->next;
			CP->next = NULL;
			ActivePolygon *AP = tActivePolygonHead.next;
			if (!AP)	tActivePolygonHead.next = (ActivePolygon*)CP;
			else{
				while (AP->next)AP = AP->next;
				AP->next = (ActivePolygon*)CP;
			}
			ClassifiedEdge *l, *r;
			findEdge(&l, y, CP->id);
			findEdge(&r, y, CP->id);
			if (!l||!r) {
				//cout << __LINE__ << ":find a polygon, but not find it's corresponding l & r edges." << endl;
				ClassifiedPolygon* next_CP = CP->next;
				CP->next = NULL;
				CP = next_CP;
				continue;
			}
			if (l->x > r->x || (l->x == r->x && l->dx > r->dx)) {
				ClassifiedEdge* tmp = l;
				l = r;
				r = tmp;
			}//使左右边位置正确
			ActiveEdgePair* AEP = (ActiveEdgePair*)malloc(sizeof(ActiveEdgePair));
			AEP->xl = l->x;
			AEP->dxl = l->dx;
			AEP->dyl = l->dy;
			AEP->xr = r->x;
			AEP->dxr = r->dx;
			AEP->dyr = r->dy;
			AEP->zl = (-CP->d - l->x * CP->a - y * CP->b) / CP->c;
			AEP->dzx = (-CP->a) / CP->c;
			AEP->dzy = CP->b / CP->c;
			AEP->id = l->id;
			AEP->color = CP->color;
			AEP->next = NULL;
			if (tActiveEdgePairHead.next==NULL) {
				tActiveEdgePairHead.next = AEP;
			}
			else{
				ActiveEdgePair* AEP_tmp = tActiveEdgePairHead.next;
				if (AEP_tmp == NULL)	tActiveEdgePairHead.next = AEP;
				else{
					while (AEP_tmp->next)	AEP_tmp = AEP_tmp->next;
					AEP_tmp->next = AEP;
				}
			}
			//画出多边形的上边界线
			GLdouble zx = AEP->zl;
			int x = AEP->xl + ((AEP->xl < 0) ? -0.5 : 0.5);
			while (x<0){
				zx += AEP->dzx;
				x++;
			}
			while (x<AEP->xr){
				if (x<screenWidth&&zx>g_zbuffer[x]) {
					g_zbuffer[x] = zx;
					Triple<GLfloat>bgColor(0.2f, 0.2f, 0.2f);
					setPixColor(x, y, bgColor);
				}
				zx += AEP->dzx;
				x++;
			}
			CP = CP_tmp;
		}
	}
}

void deepthUpdate(GLint y) {
	ActiveEdgePair* AEP = tActiveEdgePairHead.next;
	while (AEP){
		GLfloat zx = AEP->zl;
		int x = AEP->xl + ((AEP->xl < 0) ? -0.5 : 0.5);
		while (x<0){
			zx += AEP->dzx;
			x++;
		}
		if (x<screenWidth&&zx>g_zbuffer[x]) {
			g_zbuffer[x] = zx;
			Triple<GLfloat> bgColor(0.2f, 0.2f, 0.2f);
			setPixColor(x,y,bgColor);
		}
		zx += AEP->dzx;
		x++;
		while (x<screenWidth&&x<AEP->xr){
			if (zx>g_zbuffer[x]) {
				g_zbuffer[x] = zx;
				Triple<GLfloat> c(1.0f, 0.0f, 0.0f);
				setPixColor(x, y, c); 
			}
			zx += AEP->dzx;
			x++;
		}
		AEP = AEP->next;
	}
}

void activeEdgeTableUpdate(GLint y) {
	ActiveEdgePair* AEP_pre = &tActiveEdgePairHead;
	ActiveEdgePair* AEP = AEP_pre->next;
	while (AEP){
		AEP->dyl--;
		AEP->dyr--;
		if (AEP->dyl<=0&&AEP->dyr<=0) {//两边同时扫描完
			ClassifiedEdge *l, *r;
			findEdge(&l,y,AEP->id);
			if (l == NULL) {
				AEP_pre->next = AEP->next;
				free(AEP);
				AEP = AEP_pre->next;
			}else{
				findEdge(&r, y, AEP->id);
				if ( l->x > r->x || (l->x == r->x && l->dx>r->dx)) {
					ClassifiedEdge* _t = l;
					l = r;
					r = _t;
				}
				ActivePolygon *AP= tActivePolygonHead.next;
				while (AP && AP->id != AEP->id) AP = AP->next;
				AEP->xl = l->x;
				AEP->dxl = l->dx;
				AEP->dyl = l->dy;
				AEP->xr = r->x;
				AEP->dxr = r->dx;
				AEP->dyr = r->dy;
				AEP->zl = (-AP->d - l->x * AP->a - y * AP->b) / AP->c;
				AEP->dzx = (-AP->a) / AP->c;
				AEP->dzy = AP->b / AP->c;

				AEP_pre = AEP;
				AEP = AEP->next;
			}
		}else{
			if (AEP->dyl <= 0) {//左边扫描完  
				ClassifiedEdge* l;
				findEdge(&l, y, AEP->id);
				if (l) {//AEP不为空表示多边形仍在活化多边形行表中  
					AEP->xl = l->x;
					AEP->dxl = l->dx;
					AEP->dyl = l->dy;
				}
				else {
					AEP_pre->next = AEP->next;
					free(AEP);
					AEP = AEP_pre->next;
					continue;
				}
			}
			else {//左边未扫描完  
				AEP->xl += AEP->dxl;
				AEP->zl = AEP->zl + AEP->dxl * AEP->dzx + AEP->dzy;
			}

			if (AEP->dyr <= 0) {//右边扫描完  
				ClassifiedEdge* r;
				findEdge(&r, y, AEP->id);
				if (r) {
					AEP->xr = r->x;
					AEP->dxr = r->dx;
					AEP->dyr = r->dy;
				}
				else {
					AEP_pre->next = AEP->next;
					free(AEP);
					AEP = AEP_pre->next;
					continue;
				}
			}
			else {//右边未扫描完  
				AEP->xr += AEP->dxr;
			}
			AEP_pre = AEP;
			AEP = AEP->next;
		}
	}
}

void activePolygonTableUpdate() {
	ActivePolygon* AP = tActivePolygonHead.next;
	ActivePolygon* AP_pre = &tActivePolygonHead;
	while (AP) {
		AP->dy = AP->dy - 1;
		if (AP->dy <= 0) {
			AP_pre->next = AP->next;
			free(AP);
			AP = AP_pre->next;
		}
		else {
			AP_pre = AP;
			AP = AP_pre->next;
		}
	}
}

void findEdge(ClassifiedEdge** e, GLint y, GLint polygon_id) {
	//根据polygon_id寻找仍在活化多边形表中的活化边
	//e - 结果边的指针的地址
	//y - 当前扫描线位置
	//polygon_id - 多边形id
	*e = NULL;
	ActivePolygon* _pAP = (tActivePolygonHead.next);
	//检查多边形是否还在活化多边形表中
	while (_pAP && (_pAP->id != polygon_id || _pAP->dy <= 1)) //_pAP->dy <= 1 这样的多边形即将被移出活化表，所以不予考虑
		_pAP = _pAP->next;
	if (_pAP != NULL) {//多边形还在活化表中
		if (!tClassifiedEdge[y]) {
			//输出下列语句也有可能是正常现象，因为有可能扫描到多边形最底部时，多边形还没有移出
			//活化多边形表
#ifdef DEBUGGER
			cout << "function:findedge accur logic error(1)" << endl;
#endif
			return;
		}
		ClassifiedEdge* _pCE = tClassifiedEdge[y];
		while (_pCE && (_pCE->id != polygon_id || _pCE->used))
			_pCE = _pCE->next;
		if (!_pCE) {
#ifdef DEBUGGER
			cout << "function:findedge not find" << endl;
#endif
			return;
		}
		_pCE->used = true;
		*e = _pCE;
	}
}

void clearData() {
	ActivePolygon *p_AP = tActivePolygonHead.next, *q_AP;
	while (p_AP) {
		q_AP = p_AP->next;
		free(p_AP);
		p_AP = q_AP;
	}
	tActivePolygonHead.next = NULL;
	ActiveEdgePair *p_AEP = tActiveEdgePairHead.next, *q_AEP;
	while (p_AEP) {
		q_AEP = p_AEP->next;
		free(p_AEP);
		p_AEP = q_AEP;
	}
	tActiveEdgePairHead.next = NULL;

	ClassifiedPolygon *p_CP, *q_CP;
	for (GLuint i = 0; i < tClassifiedPolygon.size(); i++) {
		if (tClassifiedPolygon[i]) {
			p_CP = tClassifiedPolygon[i];
			while (p_CP) {
				q_CP = p_CP->next;
				free(p_CP);
				p_CP = q_CP;
			}
			tClassifiedPolygon[i] = NULL;
		}
	}
	ClassifiedEdge *p_CE, *q_CE;
	for (GLuint i = 0; i < tClassifiedEdge.size(); i++) {
		if (tClassifiedEdge[i]) {
			p_CE = tClassifiedEdge[i];
			while (p_CE) {
				q_CE = p_CE->next;
				free(p_CE);
				p_CE = q_CE;
			}
			tClassifiedEdge[i] = NULL;
		}
	}
}
Point roundVertex(Triple<GLfloat> v) {//计算浮点顶点v的最近邻整点p
	Point p;
	p.x = v.x + ((v.x < 0) ? -0.5 : 0.5);
	p.y = v.y + ((v.y < 0) ? -0.5 : 0.5);
	p.z = v.z;
	return p;
}

void loadModel(string path) {
	// 定义一个导入器 
	Assimp::Importer import;
	const aiScene* scene = import.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs);//转换所有的模型的原始几何形状为三角形

	if (!scene || scene->mFlags == AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
	{
		cout << "ERROR::ASSIMP::" << import.GetErrorString() << endl;
		return;
	}
	processNode(scene->mRootNode,scene);
}
void processNode(aiNode* node, const aiScene* scene)
{

	glm::mat4 model;
	glm::mat4 view = camera.GetViewMatrix();
	model = glm::translate(model, glm::vec3(300.0f, 300.0f, 0.0f));
	GLfloat angle = (GLfloat)glfwGetTime();
	model = glm::rotate(model, angle, glm::vec3(1.0f, 1.0f, 0.0f));
	glm::mat4 MVP = view * model;
	// Process all the node's meshes (if any)
	for (GLuint i = 0; i < node->mNumMeshes; i++)
	{
		aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
		for (GLint i = 0; i < mesh->mNumVertices; i++) {
			Triple<GLfloat> Point_tmp;
			glm::vec4 p(1.0f, 1.0f, 1.0f,1.0f);
			p.x = mesh->mVertices[i].x*scale;
			p.y = mesh->mVertices[i].y*scale;
			p.z = mesh->mVertices[i].z*scale;
			p = MVP*p;
			Point_tmp.x = p.x; 
			Point_tmp.y = p.y;
			Point_tmp.z = p.z;
			//cout << Point_tmp.x << " " << Point_tmp.y << " " << Point_tmp.z << endl;
			vertex.push_back(Point_tmp);
		}
		for (GLuint i = 0; i < mesh->mNumFaces; i++)
		{
			aiFace face_tmp = mesh->mFaces[i];
			vector<Triple<GLfloat>> v;
			for (GLuint j = 0; j < face_tmp.mNumIndices; j++) {
				//cout << face_tmp.mIndices[j] << endl;
				v.push_back(vertex[face_tmp.mIndices[j]]);
			}
			face.push_back(v);
			//cout << face_tmp.mNumIndices << endl;
		}
	}

	// Then do the same for each of its children
	for (GLuint i = 0; i < node->mNumChildren; i++)
	{
		processNode(node->mChildren[i], scene);
	}
}

#pragma endregion