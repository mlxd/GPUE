
/*
Unused. Future extension.
*/
namespace BEC2D{
	class Wavefunction{
		private:
			int2 gridSize;
			double2 dimMax;
			double2 *wfc = new double2[xDim*yDim];
		
		public:
			Wavefunction();
			Wavefunction(int xDim, int yDim, double xMax, double yMax);
			bool setGridSize(int xDim, int yDim);
			int2 getGridSize(int xDim, int yDim);
			double2 initWfc();
			double2 &getWfc();
	}
	
	Wavefunction::Wavefunction(){
		
	}
	Wavefunction::Wavefunction(int xDim, int yDim, double xMax, double yMax){
		
	}
	Wavefunction::setGridSize(int xDim, int yDim){
		
	}
}


