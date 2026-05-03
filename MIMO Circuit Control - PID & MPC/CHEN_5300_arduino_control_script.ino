// Analog pins use on UNO board
int analogPin1 = 3;     
int analogPin2 = 10;     
int analogPin3 = 11;
//int analogPin4 = 6;

// Variable setup
int data = 0;           
String userInput;
float values[3];
int v_out[3];
int vi[5];
char *ptr;

// Setup serial port
void setup(){

  Serial.begin(9600);
  pinMode(analogPin1,OUTPUT);
  pinMode(analogPin2,OUTPUT);
  pinMode(analogPin3,OUTPUT);
  //pinMode(analogPin4,OUTPUT);

}

// While port is active, read and write voltage data
void loop(){
if(Serial.available() > 0){ 

      // Read Python input
      userInput = Serial.readStringUntil('\n');      

      // Turn string into C character array
      char charBuf[userInput.length() + 1];
      userInput.toCharArray(charBuf, userInput.length() + 1);

      // Separate string array into floats
      ptr = strtok(charBuf, ","); // Split at the comma
      int i = 0;
      while (ptr != NULL && i < 3) {
        values[i] = atof(ptr); // Convert string token to float
        ptr = strtok(NULL, ",");
        i++;
      }
      
      // Write Python voltage inputs to respective pins
      analogWrite(analogPin1,(int)values[0]);
      analogWrite(analogPin2,(int)values[1]);
      analogWrite(analogPin3,(int)values[2]);
      //analogWrite(analogPin4,(int)values[3]);

      // Read voltage outputs
      vi[4] = analogRead(A0);
      vi[2] = analogRead(A3);
      vi[0] = analogRead(A1);
      vi[1] = analogRead(A2);
      vi[3] = analogRead(A4);

      // Voltage across selected components
      v_out[0] = vi[4] - vi[2];
      v_out[1] = vi[0];
      v_out[2] = vi[1] - vi[3];
    
      // Print voltage values to serial back to Python
      for (int i = 0; i < 3; i++) {
        Serial.print(v_out[i]*(5.0/1023.0)); // Convert bits to voltage
        Serial.print(","); // Add space between values
      }
      Serial.print("\n"); // Move to next line after printing array
  
  } // Serial.available
} // Void Loop