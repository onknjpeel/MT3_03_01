#include <Novice.h>
#include <Vector3.h>
#include <imgui.h>
#include <Matrix4x4.h>
#include <cstdint>
#include <corecrt_math_defines.h>
#include <corecrt_math.h>
#include <cmath>

const char kWindowTitle[] = "LC1B_06_オノ_ケンジピエール_タイトル";

static const int kWindowWidth = 1280;
static const int kWindowHeight = 720;

struct Sphere {
	Vector3 center;
	float radius;
};

float cot(float x) {
	float cot;
	cot = 1 / tanf(x);
	return cot;
}

const float pi = float(M_PI);

Vector3 Multiply(float scalar, const Vector3& v) {
	Vector3 result;
	result.x = v.x * scalar;
	result.y = v.y * scalar;
	result.z = v.z * scalar;
	return result;
}

Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 result;
	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;
	result.z = v1.z + v2.z;
	return result;
}

Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	Vector3 result;
	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;
	result.z = v1.z - v2.z;
	return result;
}

Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	Matrix4x4 result;
	result = {
		(1 / aspectRatio) * cot(fovY / 2),0,0,0,
		0,cot(fovY / 2),0,0,
		0,0,farClip / (farClip - nearClip),1,
		0,0,(-nearClip * farClip) / (farClip - nearClip),0
	};
	return result;
}

Matrix4x4 MakeVieportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {
	Matrix4x4 result;
	result = {
		width / 2,0,0,0,
		0,-height / 2,0,0,
		0,0,maxDepth - minDepth,0,
		left + (width / 2),top + (height / 2),minDepth,1
	};
	return result;
}

Matrix4x4 Add(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;
	result = {
		m1.m[0][0] + m2.m[0][0],m1.m[0][1] + m2.m[0][1],m1.m[0][2] + m2.m[0][2],m1.m[0][3] + m2.m[0][3],
		m1.m[1][0] + m2.m[1][0],m1.m[1][1] + m2.m[1][1],m1.m[1][2] + m2.m[1][2],m1.m[1][3] + m2.m[1][3],
		m1.m[2][0] + m2.m[2][0],m1.m[2][1] + m2.m[2][1],m1.m[2][2] + m2.m[2][2],m1.m[2][3] + m2.m[2][3],
		m1.m[3][0] + m2.m[3][0],m1.m[3][1] + m2.m[3][1],m1.m[3][2] + m2.m[3][2],m1.m[3][3] + m2.m[3][3]
	};
	return result;
}

Matrix4x4 Subtract(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;
	result = {
		m1.m[0][0] - m2.m[0][0],m1.m[0][1] - m2.m[0][1],m1.m[0][2] - m2.m[0][2],m1.m[0][3] - m2.m[0][3],
		m1.m[1][0] - m2.m[1][0],m1.m[1][1] - m2.m[1][1],m1.m[1][2] - m2.m[1][2],m1.m[1][3] - m2.m[1][3],
		m1.m[2][0] - m2.m[2][0],m1.m[2][1] - m2.m[2][1],m1.m[2][2] - m2.m[2][2],m1.m[2][3] - m2.m[2][3],
		m1.m[3][0] - m2.m[3][0],m1.m[3][1] - m2.m[3][1],m1.m[3][2] - m2.m[3][2],m1.m[3][3] - m2.m[3][3]
	};
	return result;
}

Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				result.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
	}
	return result;
}

Matrix4x4 MakeRotateMatrix(const Vector3& rotate) {
	Matrix4x4 rotateX;
	rotateX = {
	1,0,0,0,
		0,std::cos(rotate.x),std::sin(rotate.x),0,
		0,-std::sin(rotate.x),std::cos(rotate.x),0,
		0,0,0,1
	};
	Matrix4x4 rotateY;
	rotateY = {
	std::cos(rotate.y),0,-std::sin(rotate.y),0,
		0,1,0,0,
		std::sin(rotate.y),0,std::cos(rotate.y),
		0,0,0,0,1
	};
	Matrix4x4 rotateZ;
	rotateZ = {
	std::cos(rotate.z),std::sin(rotate.z),0,0,
		-std::sin(rotate.z),std::cos(rotate.z),0,0,
		0,0,1,0,
		0,0,0,1
	};
	Matrix4x4 result = Multiply(rotateX, Multiply(rotateY, rotateZ));
	return result;
}

Matrix4x4 Inverse(const Matrix4x4& m) {
	float determinant =
		+m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3]
		+ m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1]
		+ m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2]

		- m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1]
		- m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3]
		- m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2]

		- m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3]
		- m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1]
		- m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2]

		+ m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1]
		+ m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3]
		+ m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2]

		+ m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3]
		+ m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1]
		+ m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2]

		- m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1]
		- m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3]
		- m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2]

		- m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0]
		- m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0]
		- m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0]

		+ m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0]
		+ m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0]
		+ m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];

	Matrix4x4 result = {};
	float recpDeterminant = 1.0f / determinant;
	result.m[0][0] = (m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] +
		m.m[1][3] * m.m[2][1] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][1] -
		m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[1][1] * m.m[2][3] * m.m[3][2]) * recpDeterminant;
	result.m[0][1] = (-m.m[0][1] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[2][3] * m.m[3][1] -
		m.m[0][3] * m.m[2][1] * m.m[3][2] + m.m[0][3] * m.m[2][2] * m.m[3][1] +
		m.m[0][2] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[2][3] * m.m[3][2]) * recpDeterminant;
	result.m[0][2] = (m.m[0][1] * m.m[1][2] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][1] +
		m.m[0][3] * m.m[1][1] * m.m[3][2] - m.m[0][3] * m.m[1][2] * m.m[3][1] -
		m.m[0][2] * m.m[1][1] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[3][2]) * recpDeterminant;
	result.m[0][3] = (-m.m[0][1] * m.m[1][2] * m.m[2][3] - m.m[0][2] * m.m[1][3] * m.m[2][1] -
		m.m[0][3] * m.m[1][1] * m.m[2][2] + m.m[0][3] * m.m[1][2] * m.m[2][1] +
		m.m[0][2] * m.m[1][1] * m.m[2][3] + m.m[0][1] * m.m[1][3] * m.m[2][2]) * recpDeterminant;

	result.m[1][0] = (-m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[1][2] * m.m[2][3] * m.m[3][0] -
		m.m[1][3] * m.m[2][0] * m.m[3][2] + m.m[1][3] * m.m[2][2] * m.m[3][0] +
		m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[1][0] * m.m[2][3] * m.m[3][2]) * recpDeterminant;
	result.m[1][1] = (m.m[0][0] * m.m[2][2] * m.m[3][3] + m.m[0][2] * m.m[2][3] * m.m[3][0] +
		m.m[0][3] * m.m[2][0] * m.m[3][2] - m.m[0][3] * m.m[2][2] * m.m[3][0] -
		m.m[0][2] * m.m[2][0] * m.m[3][3] - m.m[0][0] * m.m[2][3] * m.m[3][2]) * recpDeterminant;
	result.m[1][2] = (-m.m[0][0] * m.m[1][2] * m.m[3][3] - m.m[0][2] * m.m[1][3] * m.m[3][0] -
		m.m[0][3] * m.m[1][0] * m.m[3][2] + m.m[0][3] * m.m[1][2] * m.m[3][0] +
		m.m[0][2] * m.m[1][0] * m.m[3][3] + m.m[0][0] * m.m[1][3] * m.m[3][2]) * recpDeterminant;
	result.m[1][3] = (m.m[0][0] * m.m[1][2] * m.m[2][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] +
		m.m[0][3] * m.m[1][0] * m.m[2][2] - m.m[0][3] * m.m[1][2] * m.m[2][0] -
		m.m[0][2] * m.m[1][0] * m.m[2][3] - m.m[0][0] * m.m[1][3] * m.m[2][2]) * recpDeterminant;

	result.m[2][0] = (m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] +
		m.m[1][3] * m.m[2][0] * m.m[3][1] - m.m[1][3] * m.m[2][1] * m.m[3][0] -
		m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][1]) * recpDeterminant;
	result.m[2][1] = (-m.m[0][0] * m.m[2][1] * m.m[3][3] - m.m[0][1] * m.m[2][3] * m.m[3][0] -
		m.m[0][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[2][1] * m.m[3][0] +
		m.m[0][1] * m.m[2][0] * m.m[3][3] + m.m[0][0] * m.m[2][3] * m.m[3][1]) * recpDeterminant;
	result.m[2][2] = (m.m[0][0] * m.m[1][1] * m.m[3][3] + m.m[0][1] * m.m[1][3] * m.m[3][0] +
		m.m[0][3] * m.m[1][0] * m.m[3][1] - m.m[0][3] * m.m[1][1] * m.m[3][0] -
		m.m[0][1] * m.m[1][0] * m.m[3][3] - m.m[0][0] * m.m[1][3] * m.m[3][1]) * recpDeterminant;
	result.m[2][3] = (-m.m[0][0] * m.m[1][1] * m.m[2][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] -
		m.m[0][3] * m.m[1][0] * m.m[2][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] +
		m.m[0][1] * m.m[1][0] * m.m[2][3] + m.m[0][0] * m.m[1][3] * m.m[2][1]) * recpDeterminant;

	result.m[3][0] = (-m.m[1][0] * m.m[2][1] * m.m[3][2] - m.m[1][1] * m.m[2][2] * m.m[3][0] -
		m.m[1][2] * m.m[2][0] * m.m[3][1] + m.m[1][2] * m.m[2][1] * m.m[3][0] +
		m.m[1][1] * m.m[2][0] * m.m[3][2] + m.m[1][0] * m.m[2][2] * m.m[3][1]) * recpDeterminant;
	result.m[3][1] = (m.m[0][0] * m.m[2][1] * m.m[3][2] + m.m[0][1] * m.m[2][2] * m.m[3][0] +
		m.m[0][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[2][1] * m.m[3][0] -
		m.m[0][1] * m.m[2][0] * m.m[3][2] - m.m[0][0] * m.m[2][2] * m.m[3][1]) * recpDeterminant;
	result.m[3][2] = (-m.m[0][0] * m.m[1][1] * m.m[3][2] - m.m[0][1] * m.m[1][2] * m.m[3][0] -
		m.m[0][2] * m.m[1][0] * m.m[3][1] + m.m[0][2] * m.m[1][1] * m.m[3][0] +
		m.m[0][1] * m.m[1][0] * m.m[3][2] + m.m[0][0] * m.m[1][2] * m.m[3][1]) * recpDeterminant;
	result.m[3][3] = (m.m[0][0] * m.m[1][1] * m.m[2][2] + m.m[0][1] * m.m[1][2] * m.m[2][0] +
		m.m[0][2] * m.m[1][0] * m.m[2][1] - m.m[0][2] * m.m[1][1] * m.m[2][0] -
		m.m[0][1] * m.m[1][0] * m.m[2][2] - m.m[0][0] * m.m[1][2] * m.m[2][1]) * recpDeterminant;

	return result;
}

Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {
	Matrix4x4 result;
	Matrix4x4 rotateM = MakeRotateMatrix(rotate);
	result = {
		scale.x * rotateM.m[0][0],scale.x * rotateM.m[0][1],scale.x * rotateM.m[0][2],0,
		scale.y * rotateM.m[1][0],scale.y * rotateM.m[1][1],scale.y * rotateM.m[1][2],0,
		scale.z * rotateM.m[2][0],scale.z * rotateM.m[2][1],scale.z * rotateM.m[2][2],0,
		translate.x,translate.y,translate.z,1
	};
	return result;
}

Vector3 Transform(Vector3 vector, Matrix4x4 matrix) {
	Vector3 result;
	result = {
		vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + matrix.m[3][0],
		vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + matrix.m[3][1],
		vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + matrix.m[3][2]
	};
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + matrix.m[3][3];

	assert(w != 0.0f);
	result.x /= w;
	result.y /= w;
	result.z /= w;

	return result;
}

#pragma region 演算子オーバーロード
Vector3 operator+(const Vector3& v1, const Vector3& v2) { return Add(v1, v2); }
Vector3 operator-(const Vector3& v1, const Vector3& v2) { return Subtract(v1, v2); }
Vector3 operator*(float s, const Vector3& v) { return Multiply(s, v); }
Vector3 operator*(const Vector3& v, float s) { return s * v; }
Vector3 operator/(const Vector3& v, float s) { return Multiply(1.0f / s, v); }
Matrix4x4 operator+(const Matrix4x4& m1, const Matrix4x4& m2) { return Add(m1, m2); }
Matrix4x4 operator-(const Matrix4x4& m1, const Matrix4x4& m2) { return Subtract(m1, m2); }
Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2) { return Multiply(m1, m2); }

Vector3 operator-(const Vector3& v) { return { -v.x,-v.y,-v.z }; }
Vector3 operator+(const Vector3& v) { return v; }

Vector3& operator*=(Vector3& v, const float& s) {
	v.x *= s;
	v.y *= s;
	v.z *= s;
	Vector3& result = v;
	return result;
}

Vector3& operator-=(Vector3& v1, const Vector3& v2) {
	v1.x -= v2.x;
	v1.y -= v2.y;
	v1.z -= v2.z;
	Vector3& result = v1;
	return result;
};

Vector3& operator+=(Vector3& v1, const Vector3& v2) {
	v1.x += v2.x;
	v1.y += v2.y;
	v1.z += v2.z;
	Vector3& result = v1;
	return result;
};

Vector3& operator/=(Vector3& v, const float& s) {
	v.x /= s;
	v.y /= s;
	v.z /= s;
	Vector3& result = v;
	return result;
};
#pragma endregion

void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {

	const float kGridHalfWidth = 2.0f;
	const uint32_t kSubdivision = 10;
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);

	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		float x = -kGridHalfWidth + (xIndex * kGridEvery);
		unsigned int color = 0xAAAAAAFF;

		Vector3 start{ x,0.0f,-kGridHalfWidth };
		Vector3 end{ x,0.0f,kGridHalfWidth };

		Vector3 startScreen = Transform(Transform(start, viewProjectionMatrix), viewportMatrix);
		Vector3 endScreen = Transform(Transform(end, viewProjectionMatrix), viewportMatrix);

		if (x == 0.0f) {
			color = BLACK;
		}
		Novice::DrawLine(int(startScreen.x), int(startScreen.y), int(endScreen.x), int(endScreen.y), color);
	}
	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {
		float z = -kGridHalfWidth + (zIndex * kGridEvery);
		unsigned int color = 0xAAAAAAFF;

		Vector3 start{ -kGridHalfWidth,0.0f,z };
		Vector3 end{ kGridHalfWidth,0.0f,z };

		Vector3 startScreen = Transform(Transform(start, viewProjectionMatrix), viewportMatrix);
		Vector3 endScreen = Transform(Transform(end, viewProjectionMatrix), viewportMatrix);

		if (z == 0.0f) {
			color = BLACK;
		}
		Novice::DrawLine(int(startScreen.x), int(startScreen.y), int(endScreen.x), int(endScreen.y), color);
	}
}

void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {

	const uint32_t kSubdivision = 20;
	const float kLonEvery = pi * 2.0f / float(kSubdivision);
	const float kLatEvery = pi / float(kSubdivision);

	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {
		float lat = -pi / 2.0f + kLatEvery * latIndex;

		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {
			float lon = lonIndex * kLonEvery;

			Vector3 a = { sphere.center.x + sphere.radius * (cosf(lat) * cosf(lon)),sphere.center.y + sphere.radius * (sinf(lat)),sphere.center.z + sphere.radius * (cosf(lat) * sinf(lon)) };
			Vector3 b = { sphere.center.x + sphere.radius * (cosf(lat + kLatEvery) * cosf(lon)),sphere.center.y + sphere.radius * (sinf(lat + kLatEvery)),sphere.center.z + sphere.radius * (cosf(lat + kLatEvery) * sinf(lon)) };
			Vector3 c = { sphere.center.x + sphere.radius * (cosf(lat) * cosf(lon + kLonEvery)),sphere.center.y + sphere.radius * (sinf(lat)),sphere.center.z + sphere.radius * (cosf(lat) * sinf(lon + kLonEvery)) };

			Vector3 worldA = Transform(Transform(a, viewProjectionMatrix), viewportMatrix);
			Vector3 worldB = Transform(Transform(b, viewProjectionMatrix), viewportMatrix);
			Vector3 worldC = Transform(Transform(c, viewProjectionMatrix), viewportMatrix);

			Novice::DrawLine(int(worldA.x), int(worldA.y), int(worldB.x), int(worldB.y), color);//LineAB
			Novice::DrawLine(int(worldA.x), int(worldA.y), int(worldC.x), int(worldC.y), color);//LineAC
		}
	}
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	Vector3 translates[3]{
		{0.2f,1.0f,0.0f},
		{0.4f,0.0f,0.0f},
		{0.3f,0.0f,0.0f}
	};

	Vector3 rotates[3]{
		{0.0f,0.0f,-6.8f},
		{0.0f,0.0f,-1.4f},
		{0.0f,0.0f,0.0f}
	};

	Vector3 scales[3]{
		{1.0f,1.0f,1.0f},
		{1.0f,1.0f,1.0f},
		{1.0f,1.0f,1.0f}
	};

	Vector3 rotate{ 0.0f,0.0f,0.0f };
	Vector3 translate{ 0.0f,0.0f,0.0f };

	Sphere sphere[3];

	Vector3 cameraTranslate{ 0.0f,1.9f,-6.49f };
	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };

	Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, rotate, translate);
	Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);
	Matrix4x4 viewMatrix = Inverse(cameraMatrix);
	Matrix4x4 viewProjectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
	Matrix4x4 worldViewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, viewProjectionMatrix));
	Matrix4x4 viewportMatrix = MakeVieportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

	Matrix4x4 shoulderWM;
	Matrix4x4 elbowWM;
	Matrix4x4 handWM;

	Matrix4x4 shoulderWVP;
	Matrix4x4 elbowWVP;
	Matrix4x4 handWVP;

	Vector3 shoulderPoint;
	Vector3 elbowPoint;
	Vector3 handPoint;

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		for (int32_t i = 0; i < 3; i++) {
			sphere[i] = { translate,0.05f};
		}

		shoulderWM = MakeAffineMatrix(scales[0], rotates[0], translates[0]);
		elbowWM = MakeAffineMatrix(scales[1], rotates[1], translates[1]);
		handWM = MakeAffineMatrix(scales[2], rotates[2], translates[2]);

		shoulderWVP = Multiply(shoulderWM, Multiply(viewMatrix, viewProjectionMatrix));
		elbowWVP = Multiply(Multiply(elbowWM,shoulderWM), Multiply(viewMatrix, viewProjectionMatrix));
		handWVP = Multiply(Multiply(handWM,Multiply(elbowWM,shoulderWM)), Multiply(viewMatrix, viewProjectionMatrix));

		shoulderPoint = Transform(Transform(sphere[0].center, shoulderWVP), viewportMatrix);
		elbowPoint = Transform(Transform(sphere[1].center, elbowWVP), viewportMatrix);
		handPoint = Transform(Transform(sphere[2].center, handWVP), viewportMatrix);

		ImGui::Begin("Windows");
		ImGui::DragFloat3("shoulder.translate", &translates[0].x, 0.01f);
		ImGui::DragFloat3("shoulder.rotate", &rotates[0].x, 0.01f);
		ImGui::DragFloat3("shoulder.scale", &scales[0].x, 0.01f);
		ImGui::DragFloat3("elbow.translate", &translates[1].x, 0.01f);
		ImGui::DragFloat3("elbow.rotate", &rotates[1].x, 0.01f);
		ImGui::DragFloat3("elbow.scale", &scales[1].x, 0.01f);
		ImGui::DragFloat3("hand.translate", &translates[2].x, 0.01f);
		ImGui::DragFloat3("hand.rotate", &rotates[2].x, 0.01f);
		ImGui::DragFloat3("hand.scale", &scales[2].x, 0.01f);
		ImGui::End();

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		DrawGrid(worldViewProjectionMatrix, viewportMatrix);

		Novice::DrawLine(int(shoulderPoint.x), int(shoulderPoint.y), int(elbowPoint.x), int(elbowPoint.y), WHITE);
		Novice::DrawLine(int(elbowPoint.x), int(elbowPoint.y), int(handPoint.x), int(handPoint.y), WHITE);

		DrawSphere(sphere[0], shoulderWVP, viewportMatrix, RED);
		DrawSphere(sphere[1], elbowWVP, viewportMatrix, GREEN);
		DrawSphere(sphere[2], handWVP, viewportMatrix, BLUE);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
