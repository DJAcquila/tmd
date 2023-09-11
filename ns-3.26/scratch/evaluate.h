typedef struct _results {
	int m_lostPackets;		//Total packets lost
	int m_rxPackets;		//Total packets received
	double m_throughput;	//Average throughput
	double m_delayAvg;		//Average end-to-end time
	double m_jitterAvg;		//Average end-to-end jitter
	int m_MacTxDrop; 		//A packet has been dropped in the MAC layer before being queued for transmission.
	int m_MacRxDrop; 		//A packet has been dropped in the MAC layer after it has been passed up from the physical layer.
	int m_PhyTxDrop; 		//Trace source indicating a packet has been dropped by the device during transmission.
	int m_PhyRxDrop; 		//Trace source indicating a packet has been dropped by the device during reception.
	int m_ArpDrop;			//Trace source indicating a packet has been dropped by ARP protocol.
	int m_MacQueueDrop;
	unsigned int m_QuantityHops; //Trace source indicating a quantity hops for packet received.
	int recalculates;
	int qtdLinkUpTrue;
	int qtdLinkUpFalse;
	int totalLowQuality;
} RESULTS;
